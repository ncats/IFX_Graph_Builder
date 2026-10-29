# POUNCE sheet validation API

`POST /pounce/api/validate` accepts parsed workbook sheets as JSON. The request
contract is [pounce_validation_request.schema.json](pounce_validation_request.schema.json).
The endpoint is stateless: it validates the supplied sheets and returns a report;
it does not create an upload or email submission session. The HTML upload and
JSON endpoint use the same workbook parser, content validators, and issue
formatting logic.

## Request

```json
{
  "schema_version": 1,
  "project": {
    "source_id": "project.xlsx",
    "sheets": [
      {
        "name": "ProjectMeta",
        "columns": ["NCATSDPI_Variable_Name", "Submitter_Value"],
        "rows": [["project_name", "Example project"]]
      }
    ]
  },
  "experiments": [],
  "stats_results": []
}
```

This short example shows the shape, not a complete valid submission. Supply
every sheet, its exact submitted name and column names, and all rows, including
the data matrices. Each row is an array with one cell per column, in column
order. Empty sheets still have `columns` and `rows: []`. Use JSON strings,
numbers, booleans, or `null` for cells; use ISO 8601 strings for dates. Blank
Excel cells may be sent as `""` or `null`. Keep repeated map rows rather than
collapsing them into an object.

`source_id` is a unique, caller-chosen workbook label without path separators
within the request. It appears in validation issues. Pair `experiments[i]` with
`stats_results[i]`; the arrays must have the same length. Use `null` for an
experiment without StatsResults, including a missing slot between two present
workbooks. Do not remove empty slots or later stats files will pair with the
wrong experiment. The server detects v2 or v3 separately
for each workbook from its map-sheet headers and normalizes both into the
shared parser. Do not add a template-version field to the request.

## Response

HTTP 200 means the request was processed. Validation problems appear in
`issues` and are counted in `error_count` and `warning_count`:

```json
{
  "schema_version": 1,
  "files": [
    {"role": "project", "source_id": "project.xlsx", "template_version": 3}
  ],
  "error_count": 1,
  "warning_count": 0,
  "issues": [
    {
      "severity": "error",
      "entity": "project",
      "field": "project_name",
      "message": "Project name is required",
      "sheet": "ProjectMeta",
      "column": "project_name",
      "row": null,
      "source_id": "project.xlsx"
    }
  ],
  "coverage": []
}
```

The issue above illustrates response fields; its message is not a promised
literal validator message. `source_id`, `sheet`, `column`, and `row` can be
`null` if an issue has no precise location. `coverage` contains analyte type,
sheet, total, mapped, unmapped IDs, and mapped percentage when a resolver is
configured. An empty coverage array does not mean that all IDs mapped.

HTTP 422 means the request or workbook format cannot be processed, such as
duplicate source IDs, unequal experiment/stats counts, malformed row lengths,
or ambiguous template version. The response uses FastAPI's `detail` field.
Requests larger than 50 MiB return HTTP 413. The response's `files` list
contains only supplied workbooks, with no entry for a `null` stats slot.

## Scope

The API reports the same validation findings as the workbook upload, including
sheet structure, mapped columns, matrix references, and content rules. It does
not run the graph ETL or submit files by email.
