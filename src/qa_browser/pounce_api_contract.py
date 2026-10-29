"""Public JSON contract for validating already decoded POUNCE workbook sheets."""

from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator


Cell = str | int | float | bool | None
Label = Annotated[str, Field(min_length=1)]
SourceId = Annotated[str, Field(min_length=1, pattern=r"^[^/\\]+$")]


class SheetTable(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)

    name: Label
    columns: list[Label]
    rows: list[list[Cell]]

    @model_validator(mode="after")
    def check_shape(self):
        if len(self.columns) != len({column.strip() for column in self.columns}):
            raise ValueError(f"{self.name}: duplicate column names")
        if any(len(row) != len(self.columns) for row in self.rows):
            raise ValueError(f"{self.name}: every row must have {len(self.columns)} cells")
        return self


class WorkbookTables(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)

    source_id: SourceId
    sheets: list[SheetTable]

    @model_validator(mode="after")
    def check_sheet_names(self):
        names = [sheet.name for sheet in self.sheets]
        if len(names) != len(set(names)):
            raise ValueError(f"{self.source_id}: duplicate sheet names")
        return self


class PounceValidationRequest(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)

    schema_version: Literal[1]
    project: WorkbookTables
    experiments: list[WorkbookTables]
    stats_results: list[WorkbookTables | None]

    @model_validator(mode="after")
    def check_workbooks(self):
        workbooks = [self.project, *self.experiments, *(item for item in self.stats_results if item is not None)]
        source_ids = [workbook.source_id for workbook in workbooks]
        if len(source_ids) != len(set(source_ids)):
            raise ValueError("source_id must be unique across the submission")
        if len(self.experiments) != len(self.stats_results):
            raise ValueError("experiments and stats_results must have matching lengths")
        return self
