#!/usr/bin/env Rscript

# Test script for segment_traced R wrapper function

# Source the function
source("R/segment_traced.R")

cat("=== Testing segment_traced R wrapper ===\n\n")

# Test 1: Single sequence without trace
cat("Test 1: Simple AT repeat (no trace)\n")
cat("--------------------------------------\n")
result1 <- segment_traced("ATATAT", "ATATGG", trace = FALSE)
print(result1)
cat("\n")

# Test 2: Single sequence with trace
cat("Test 2: CG repeat with trace output\n")
cat("--------------------------------------\n")
result2 <- segment_traced("CGCGCG", "CGCGAA", trace = TRUE)
cat("\nParsed result:\n")
print(result2)
cat("\n")

# Test 3: Multiple sequences
cat("Test 3: Multiple sequences (no trace)\n")
cat("--------------------------------------\n")
result3 <- segment_traced(
  string = c("ATATAT", "CGCGCG", "TTTTT", "ACGTACGT"),
  context = c("ATATGG", "CGCGAA", "TTTGGG", "ACGTACGTNN"),
  trace = FALSE
)
cat("\nResult structure:\n")
str(result3)
cat("\n")

# Test 4: Compare with original segment function (if available)
cat("Test 4: Comparison with original segment() function\n")
cat("--------------------------------------\n")

# Check if original segment is available
if (file.exists("src/code.cpp") && requireNamespace("Rcpp", quietly = TRUE)) {
  # Try to load the package
  tryCatch({
    devtools::load_all(quiet = TRUE)

    # Run original segment
    original_result <- segment(
      string = c("ATATAT", "CGCGCG"),
      context = c("ATATGG", "CGCGAA")
    )

    # Run traced segment
    traced_result <- segment_traced(
      string = c("ATATAT", "CGCGCG"),
      context = c("ATATGG", "CGCGAA"),
      trace = FALSE
    )

    # Convert traced_result list to data frame for comparison
    traced_df <- as.data.frame(traced_result, stringsAsFactors = FALSE)

    cat("Original segment() result:\n")
    print(original_result)
    cat("\n")

    cat("segment_traced() result:\n")
    print(traced_df)
    cat("\n")

    # Check if they match
    cat("Results match: ")
    all_match <- all.equal(original_result, traced_df, check.attributes = FALSE)
    if (isTRUE(all_match)) {
      cat("YES ✓\n")
    } else {
      cat("NO\n")
      cat("Differences:\n")
      print(all_match)
    }

  }, error = function(e) {
    cat("Could not load package for comparison:", e$message, "\n")
  })
} else {
  cat("Skipping comparison (package not available or Rcpp not installed)\n")
}

cat("\n=== Tests completed ===\n")
