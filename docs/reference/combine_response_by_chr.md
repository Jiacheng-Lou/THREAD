# Combine chromosome-level response files

Combine chromosome-level response files

## Usage

``` r
combine_response_by_chr(
  output_dir,
  pattern = "^chr.*_response\\.rds$",
  output_file = NULL
)
```

## Arguments

- output_dir:

  Directory containing chr\*\_response.rds files.

- pattern:

  File pattern for chromosome-level response files.

- output_file:

  Optional merged output RDS path.

## Value

A data.frame containing the merged response table.
