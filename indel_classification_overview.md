# indel_classification.R Overview

This document provides a detailed breakdown of the `R/indel_classification.R` file (323 lines), which contains the core logic for the **89-channel classification system**.

## Main Functions

### 1. `indel_classifier89(indels, genome.v, verbose=TRUE)` (lines 8-31)

The main entry point for 89-channel classification:

**Parameters:**
- `indels`: A dataframe with columns: "Sample", "chr", "position", "REF", "ALT"
- `genome.v`: Genome version ("hg19" or "hg38")
- `verbose`: Control output messages (default: TRUE)

**Process:**
1. Prepares indels using `prepare_indels()` - extracts flanking sequences
2. Segments them using `segment_indels()` - identifies repeats and microhomology
3. Calculates microhomology lengths: `mh_length = indel.length - spacer_length`
4. Caps `mh_length` at `prime3_rep_length` when it exceeds it
5. Assigns channels using `assign_channels_m5()`

**Returns:** Classified indel list with channel assignments

---

### 2. `assign_channels_m5(indel.df)` (lines 39-321)

A complex 280+ line function that assigns each indel to one of 89 channels. This is the heart of the classification system.

**Hierarchical Type System:**
- `type_1`: Basic type (Ins/Del/Complex)
- `type_2`: Color grouping (e.g., `[InsC]`, `[InsT]`, `Ins_nMer`, `Del_Spaced`)
- `type_3`: 21-channel level with repeat level (e.g., `[InsC]ShortRep_leq4`)
- `type_4`: **89-channel level** (most granular, used for final catalogues)

---

## Detailed Classification Rules

### Single Base Insertions: [+C] and [+T] (lines 81-137)

#### [+C] Insertions (lines 82-112)

**Special Case:**
- `A[Ins(C):R0]A` or `A[Ins(C):R0]T` (lines 97-101)
  - No repeats (R=0)
  - 5' flanking base is A
  - 3' flanking base is A or T

**General Cases (all other [+C]):**
- `Ins(C):R(0,3)` - 0-3 repeats (line 108)
- `Ins(C):R(4,6)` - 4-6 repeats (line 109)
- `Ins(C):R(7,9)` - 7-9 repeats (line 110)

#### [+T] Insertions (lines 114-135)

All [+T] insertions include flanking base context (5' and 3'):

- `X[Ins(T):R(0,4)]Y` - 0-4 repeats (lines 115-120)
  - X, Y can be A, C, or G
  - Creates 9 channels (3×3 combinations)

- `X[Ins(T):R(5,7)]Y` - 5-7 repeats (lines 122-128)
  - 9 channels with flanking context

- `X[Ins(T):R(8,9)]Y` - 8-9 repeats (lines 129-135)
  - 9 channels with flanking context

**Total [+T] channels:** 27 (9 + 9 + 9)

---

### Multi-base Insertions (lines 138-188)

#### N-mer with Tandem Repeats ≥2 (lines 139-148)
- **Condition:** `indel.length > 1` AND `prime3_reps >= 2` AND `spacer_length == 0`
- **Channels:**
  - `Ins(2,):R(2,4)` - 2-4 repeats
  - `Ins(2,):R(5,9)` - 5-9 repeats

#### Non-repetitive, TR=1 (lines 151-160)
- **Condition:** `indel.length > 1` AND `prime3_reps == 1` AND `spacer_length == 0`
- **Channels:**
  - `Ins(2,4):R1` - Length 2-4 bp
  - `Ins(5,):R1` - Length ≥5 bp

#### No Repeats, TR=0 (lines 164-188)
- **Condition:** `indel.length > 1` AND `prime3_reps == 0`
- Two sub-cases:
  1. `spacer_length == 0` (lines 164-173): Whole indel is spacer
  2. `spacer_length > 0` (lines 177-188): Has spacer content
- **Channels:**
  - `Ins(2,4):R0` - Length 2-4 bp
  - `Ins(5,):R0` - Length ≥5 bp

---

### Single Base Deletions: [-C] and [-T] (lines 192-262)

#### [-C] Deletions (lines 212-238)

**Specific repeat counts with A/T flanking:**
- `[Del(C):R1]A` or `[Del(C):R1]T` - 1 repeat (lines 212-217)
- `[Del(C):R2]A` or `[Del(C):R2]T` - 2 repeats
- `[Del(C):R3]A` or `[Del(C):R3]T` - 3 repeats
- `[Del(C):R(4,5)]A` or `[Del(C):R(4,5)]T` - 4-5 repeats (lines 219-225)

**Special case with G flanking:**
- `[Del(C):R(1,5)]G` - 1-5 repeats, 3' flanked by G (lines 226-232)

**High repeat count:**
- `Del(C):R(6,9)` - 6-9 repeats (lines 234-239)

**Total [-C] channels:** 9

#### [-T] Deletions (lines 241-261)

All [-T] deletions include flanking base context:

- `X[Del(T):R(1,4)]Y` - 1-4 repeats (lines 242-247)
  - 9 channels (3×3 flanking combinations)

- `X[Del(T):R(5,7)]Y` - 5-7 repeats (lines 248-253)
  - 9 channels

- `X[Del(T):R(8,9)]Y` - 8-9 repeats (lines 255-261)
  - 9 channels

**Total [-T] channels:** 27

---

### Multi-base Deletions (lines 264-305)

#### N-mer Deletions (lines 265-276)
- **Condition:** `indel.length > 1` AND `prime3_reps > 0` AND `spacer_length == 0`
- **Channels:**
  - `Del(2,8):U(1,2):R(2,4)` - Unit length 1-2, 2-4 repeats
  - `Del(2,):U(1,2):R(5,9)` - Unit length 1-2, 5-9 repeats
  - `Del(3,):U(3,):R2` - Unit length ≥3, 2 repeats
  - `Del(3,):U(3,):R(3,9)` - Unit length ≥3, 3-9 repeats

**Note:** U = unit length, R = repeat count

#### Spaced Deletions with Microhomology (lines 278-293)
- **Condition:** `indel.length > 1` AND `prime3_reps > 0` AND `spacer_length > 0`
- **Type classification:**
  - `Del_Spaced_short_leq5` - Length ≤5 bp
  - `Del_Spaced_long_g5` - Length >5 bp

**Channels:**
- Short deletions (≤5 bp):
  - `Del(2,5):M1` - Microhomology = 1
  - `Del(3,5):M2` - Microhomology = 2
  - `Del(4,5):M(3,4)` - Microhomology = 3-4

- Long deletions (≥6 bp):
  - `Del(6,):M1` - Microhomology = 1
  - `Del(6,):M2` - Microhomology = 2
  - `Del(6,):M3` - Microhomology = 3
  - `Del(6,):M(4,)` - Microhomology ≥4

**Note:** M = microhomology length

#### Non-repetitive Deletions (lines 296-305)
- **Condition:** `indel.length > 1` AND `prime3_reps == 0`
- **Channels:**
  - `Del(2,4):R1` - Length 2-4 bp
  - `Del(5,):R1` - Length ≥5 bp

---

### Complex Indels (lines 307-314)

**Condition:** `indel.type == "DI"` (simultaneous deletion and insertion)

**Channel:** All assigned to single `Complex` channel

**Processing:**
- Sets `indel.length` to length of REF
- Sets `change` to the REF sequence

---

## Key Implementation Details

### Pyrimidine Normalization
The function uses pyrimidine-normalized versions of sequences for strand-neutral classification:
- `change_pyr`: Pyrimidine version of the indel sequence
- `slice5_pyr`: Pyrimidine version of 5' flanking sequence
- `slice3_pyr`: Pyrimidine version of 3' flanking sequence

### Flanking Base Extraction
For single-base indels, the function extracts flanking bases:
- `rep_slice5`: Base immediately 5' to the repeat region
- `rep_slice3`: Base immediately 3' to the repeat region

**Different logic for forward vs reverse strand:**
- Forward (C/T): Uses end of slice5, position after repeats in slice3
- Reverse (G/A): Adjusts positions using `str_sub` negative indexing

### Repeat Level Classification
Three repeat levels used in type_3:
- `NonRep`: Insertions with <2 repeats, Deletions with 1 repeat
- `ShortRep_leq4`: 2-4 repeats (Ins) or 2-4 repeats (Del)
- `LongRep_g4`: >4 repeats

### Defensive Programming
The function checks dataframe dimensions before processing each category:
```r
if(dim(indel.b1)[1] > 0) {
  # ... processing logic
}
```
This prevents errors when a particular indel category is empty.

---

## Channel Count Summary

**Single base insertions:**
- [+C]: 4 channels (1 special + 3 general)
- [+T]: 27 channels (3 repeat groups × 9 flanking combinations)

**Multi-base insertions:**
- N-mer repeats: 2 channels
- Non-repetitive: 4 channels (2 R0 + 2 R1)

**Single base deletions:**
- [-C]: 9 channels
- [-T]: 27 channels

**Multi-base deletions:**
- N-mer: 4 channels
- Spaced: 7 channels
- Non-repetitive: 2 channels

**Complex:** 1 channel

**Total:** 89 channels

---

## Usage Example

```r
# Load indel data
indels <- data.frame(
  Sample = "sample1",
  chr = "chr1",
  position = c(100, 200, 300),
  REF = c("TC", "A", "GGCC"),
  ALT = c("T", "AT", "G")
)

# Classify indels
classified <- indel_classifier89(indels, genome.v = "hg38")

# Result includes type_4 column with 89-channel assignments
table(classified$type_4)
```

---

## Related Functions

- `prepare_indels()`: In `R/auxiliary_functions.R` - extracts flanking sequences
- `segment_indels()`: Uses C++ code in `src/code.cpp` - finds repeats
- `gen_catalogue89()`: In `R/gen_indelcatalogue.R` - generates counts per channel
- `indel_highspecific()`: In `R/indel_highspecific.R` - filters high-repeat regions
