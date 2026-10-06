"""Atomic SQLite artifact creation with exact-row deduplication and validation."""
from collections import defaultdict
import json
import os
from pathlib import Path
import re
import sqlite3
import tempfile


class SQLiteWriter:
    @staticmethod
    def validate_output(output, *, overwrite=False):
        output = Path(output)
        if output.is_symlink() or (output.exists() and not output.is_file()):
            raise ValueError(f'Output must be an ordinary file path: {output}')
        if output.exists() and not overwrite:
            raise FileExistsError(f'Output already exists: {output}. Use --overwrite to replace it, or choose another path.')
        if overwrite and any(Path(str(output) + suffix).exists() for suffix in ('-wal', '-shm', '-journal')):
            raise ValueError(f'Output has SQLite journal/WAL files: {output}. Close its database connections before replacing it.')

    def __init__(self, output, *, overwrite=False):
        self.output = Path(output)
        self.overwrite = overwrite
        self.validate_output(self.output, overwrite=overwrite)
        self.output.parent.mkdir(parents=True, exist_ok=True)
        fd, temporary = tempfile.mkstemp(prefix=f".{self.output.name}.", suffix=".tmp", dir=self.output.parent)
        os.close(fd)
        self.temporary = Path(temporary)
        self.db = sqlite3.connect(self.temporary)
        self.buffers = defaultdict(list)
        self.columns = {}
        self.indexes = []

    def initialize(self):
        statement = ""
        for line in Path(__file__).with_name("sqlite_schema.sql").read_text().splitlines(True):
            statement += line
            if sqlite3.complete_statement(statement):
                if statement.lstrip().upper().startswith("CREATE INDEX") or statement.lstrip().upper().startswith("CREATE UNIQUE INDEX"):
                    self.indexes.append(statement)
                else:
                    self.db.executescript(statement)
                statement = ""
        for (table,) in self.db.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall():
            self.columns[table] = [r[1] for r in self.db.execute(f'PRAGMA table_info("{table}")')]
            self.db.execute(f'CREATE TABLE "_raw_{table}" AS SELECT * FROM "{table}" WHERE 0')

    def add(self, table, **row):
        columns = self.columns[table]
        if set(row) != set(columns):
            raise ValueError(f"Incorrect {table} columns: missing {set(columns)-set(row)}, extra {set(row)-set(columns)}")
        self.buffers[table].append(tuple(row[c] for c in columns))
        if len(self.buffers[table]) >= 5000:
            self.flush(table)

    def flush(self, table):
        rows = self.buffers[table]
        if rows:
            marks = ",".join("?" for _ in self.columns[table])
            self.db.executemany(f'INSERT INTO "_raw_{table}" VALUES ({marks})', rows)
            rows.clear()

    def retain_tables(self, names):
        """Keep only the tables needed for a diagnostic export."""
        names = set(names)
        if not names <= self.columns.keys():
            raise ValueError(f"Unknown retained tables: {names - self.columns.keys()}")
        for table in list(self.columns):
            if table not in names:
                self.buffers.pop(table, None)
                self.db.execute(f'DROP TABLE "_raw_{table}"')
                self.db.execute(f'DROP TABLE "{table}"')
                del self.columns[table]
        self.indexes = [statement for statement in self.indexes
                        if (match := re.search(r'\bON\s+"?([A-Za-z_]\w*)"?\s*\(', statement, re.I))
                        and match.group(1) in names]

    def finish(self, metadata, verify_source):
        for table in self.columns:
            self.flush(table)
            # Exact duplicates collapse. Different rows sharing a legacy primary
            # key are a mapping error, never silently discarded by INSERT IGNORE.
            self.db.execute(f'INSERT INTO "{table}" SELECT DISTINCT * FROM "_raw_{table}"')
            self.db.execute(f'DROP TABLE "_raw_{table}"')
        for statement in self.indexes:
            self.db.executescript(statement)
        counts = {table: self.db.execute(f'SELECT count(*) FROM "{table}"').fetchone()[0]
                  for table in self.columns if table != "ramp_export_metadata"}
        metadata = {**metadata, "row_counts": counts}
        if 'ramp_export_metadata' in self.columns:
            self.db.execute("INSERT INTO ramp_export_metadata VALUES (?,?)", ("manifest", json.dumps(metadata, sort_keys=True)))
        errors = self.db.execute("PRAGMA foreign_key_check").fetchmany(10)
        if errors:
            raise ValueError(f"SQLite relationship integrity failed: {errors}")
        result = self.db.execute("PRAGMA integrity_check").fetchall()
        if result != [("ok",)]:
            raise ValueError(f"SQLite integrity check failed: {result[:10]}")
        self.db.commit()
        # Raw staging pages are removed before publication, keeping the artifact compact.
        self.db.execute("VACUUM")
        verify_source()
        self.db.close()
        with self.temporary.open("rb") as handle:
            os.fsync(handle.fileno())
        self.validate_output(self.output, overwrite=self.overwrite)
        if self.overwrite:
            os.replace(self.temporary, self.output)
        else:
            # Atomic no-clobber publication also handles concurrent exports.
            os.link(self.temporary, self.output)
            self.temporary.unlink()
        return metadata

    def close(self):
        self.db.close()
        self.temporary.unlink(missing_ok=True)
