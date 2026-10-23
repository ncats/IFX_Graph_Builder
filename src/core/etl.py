import datetime
import time
import uuid
from dataclasses import dataclass, field
from typing import List, Dict
import humanize

from src.interfaces.id_resolver import IdResolver
from src.interfaces.input_adapter import InputAdapter
from src.interfaces.output_adapter import OutputAdapter


@dataclass
class ETL:
    input_adapters: List[InputAdapter]
    output_adapters: List[OutputAdapter]
    resolver_map: Dict[str, IdResolver] = field(default_factory=dict)
    curation_snapshots: dict = field(default_factory=dict)
    post_adapters: List[InputAdapter] = field(default_factory=list)

    def create_or_truncate_datastores(self, truncate_tables: bool = None):
        for output_adapter in self.output_adapters:
            if not output_adapter.create_or_truncate_datastore(truncate_tables=truncate_tables):
                raise Exception("operation cancelled")

    def do_etl(self, do_post_processing = True, clean_edges: bool = True, resume: bool = False,
               run_id: str | None = None):
        if resume and self.curation_snapshots:
            raise RuntimeError(
                "Resume is not safe when curations are configured because a changed or removed "
                "property decision may require restoring a freshly derived baseline. Run a clean build."
            )
        total_start_time = time.time()
        effective_run_id = run_id or f"etl_{datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%SZ')}_{uuid.uuid4().hex[:8]}"
        print(f"ETL run id: {effective_run_id}")

        completed_adapters = set()
        if effective_run_id:
            if resume:
                for output_adapter in self.output_adapters:
                    completed_adapters |= output_adapter.get_completed_adapter_names(effective_run_id)
                if completed_adapters:
                    print(f"Resume mode: skipping {len(completed_adapters)} completed adapters")
            else:
                for output_adapter in self.output_adapters:
                    output_adapter.reset_run_state(effective_run_id)

        if resume and any(
            name.startswith("post:") for name in completed_adapters
        ):
            incomplete_inputs = [
                adapter.get_name()
                for adapter in self.input_adapters
                if f"input:{adapter.get_name()}" not in completed_adapters
                and adapter.get_name() not in completed_adapters
            ]
            if incomplete_inputs:
                raise RuntimeError(
                    "Cannot resume completed post adapters while primary adapters are incomplete: "
                    + ", ".join(incomplete_inputs)
                    + ". Run a fresh build."
                )

        for output_adapter in self.output_adapters:
            output_adapter.do_pre_processing()

        def run_adapter_phase(adapters: List[InputAdapter], phase: str, allow_resume_skip: bool) -> None:
            adapter_total = len(adapters)
            for adapter_position, input_adapter in enumerate(adapters, start=1):
                display_name = input_adapter.get_name()
                checkpoint_name = f"{phase}:{display_name}"
                completed = checkpoint_name in completed_adapters or (
                    phase == "input" and display_name in completed_adapters
                )
                if resume and allow_resume_skip and completed:
                    print(
                        f"Skipping completed {phase} adapter "
                        f"[{adapter_position}/{adapter_total}]: {display_name}"
                    )
                    continue

                start_time = time.time()
                print(f"Running {phase} [{adapter_position}/{adapter_total}]: {display_name}")
                for output_adapter in self.output_adapters:
                    output_adapter.mark_adapter_running(
                        run_id=effective_run_id,
                        adapter_name=checkpoint_name,
                        adapter_position=adapter_position,
                        adapter_total=adapter_total,
                    )
                count = 0
                try:
                    for resolved_list in input_adapter.get_resolved_and_provenanced_list(
                        resolver_map=self.resolver_map
                    ):
                        count += len(resolved_list)
                        for output_adapter in self.output_adapters:
                            resolved_list = output_adapter.preprocess_objects(resolved_list)
                            output_adapter.store(
                                resolved_list,
                                single_source=input_adapter.is_single_source(),
                                field_conflict_behavior=input_adapter.get_field_conflict_behavior(),
                            )
                except Exception as exc:
                    for output_adapter in self.output_adapters:
                        output_adapter.mark_adapter_failed(
                            run_id=effective_run_id,
                            adapter_name=checkpoint_name,
                            error_message=str(exc),
                            adapter_position=adapter_position,
                            adapter_total=adapter_total,
                        )
                    raise

                elapsed_time = time.time() - start_time
                print(f"\tElapsed time: {elapsed_time:.4f} seconds merging {count} records")
                for output_adapter in self.output_adapters:
                    output_adapter.flush_incremental_metadata()
                    output_adapter.mark_adapter_completed(
                        run_id=effective_run_id,
                        adapter_name=checkpoint_name,
                        records_written=count,
                        adapter_position=adapter_position,
                        adapter_total=adapter_total,
                    )

        run_adapter_phase(self.input_adapters, "input", allow_resume_skip=True)
        if self.curation_snapshots:
            for output_adapter in self.output_adapters:
                output_adapter.apply_curation_snapshots(
                    self.curation_snapshots,
                    phase="pre_post",
                )
        run_adapter_phase(self.post_adapters, "post", allow_resume_skip=True)

        if self.curation_snapshots:
            for output_adapter in self.output_adapters:
                output_adapter.apply_curation_snapshots(
                    self.curation_snapshots,
                    phase="final",
                )

        if do_post_processing:
            for output_adapter in self.output_adapters:
                output_adapter.do_post_processing(clean_edges=clean_edges)

        total_elapsed_time = time.time() - total_start_time
        elapsed_timedelta = datetime.timedelta(seconds=total_elapsed_time)

        # Format the elapsed time using humanize
        formatted_time = humanize.precisedelta(elapsed_timedelta, format='%0.0f')

        print(f"\tTotal elapsed time: {formatted_time}")
