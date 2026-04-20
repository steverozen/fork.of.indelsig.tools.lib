#!/usr/bin/env Rscript

# Test script for seg_traced (single sequence function)

# Source the function
source("R/segment_traced.R")

cat("=== Testing seg_traced (single sequence) ===\n\n")

# Test 1: Single sequence without trace
cat("Test 1: Simple AT repeat (no trace)\n")
cat("--------------------------------------\n")
result1 <- seg_traced("ATATAT", "ATATGG", trace = FALSE)
print(result1)
cat("\n")

# Test 2: Single sequence with trace
cat("Test 2: CG repeat with trace output\n")
cat("--------------------------------------\n")
result2 <- seg_traced("CGCGCG", "CGCGAA", trace = TRUE)
cat("\nParsed result:\n")
print(result2)
cat("\n")

# Test 3: Verify it's a scalar output (not vectors)
cat("Test 3: Verify output types\n")
cat("--------------------------------------\n")
cat("unit is character scalar:", is.character(result1$unit), "length:", length(result1$unit), "\n")
cat("unit_length is integer scalar:", is.integer(result1$unit_length), "length:", length(result1$unit_length), "\n")
cat("\n")

# Test 4: Verify segment_traced (vectorized) still works
cat("Test 4: Verify segment_traced (vectorized) still works\n")
cat("--------------------------------------\n")
result4 <- segment_traced(
  string = c("ATATAT", "CGCGCG"),
  context = c("ATATGG", "CGCGAA"),
  trace = FALSE
)
cat("Result structure:\n")
str(result4)
cat("\n")

cat("=== Tests completed ===\n")
