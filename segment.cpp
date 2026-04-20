#include <Rcpp.h>
using namespace std;
using namespace Rcpp;

/*
 * Comparison function for sorting segmentations by priority
 *
 * Segmentation vector format: [unit_length, internal_reps, spacer_length, prime3_reps]
 *
 * Priority order (best segmentation has):
 * 1. HIGHEST 3' context repeats [3] - most important
 * 2. HIGHEST internal repeats [1]
 * 3. LOWEST spacer length [2] - prefer clean repeats
 * 4. LOWEST unit length [0] - prefer simpler repeat units
 */
bool compareSegmentations(std::vector<int> v1, std::vector<int> v2){
    // Priority 1: HIGHEST 3' context repeats (most important)
    // Indels are more likely in regions that continue repeating into flanking sequence
    if(v1[3] > v2[3]){
        return true;   // v1 is better
    }else if(v1[3] < v2[3]){
        return false;  // v2 is better
    } else{
        // Priority 2: HIGHEST internal repeats (when 3' repeats are tied)
        // More internal repetition suggests stronger repeat-mediated mechanism
        if(v1[1] > v2[1]){
            return true;
        }else if(v1[1] < v2[1]){
            return false;
        }else{
            // Priority 3: LOWEST spacer length (when internal repeats are tied)
            // Cleaner repeat patterns (less leftover bases) are preferred
            if(v1[2] < v2[2]){
                return true;
            }else if(v1[2] > v2[2]){
                return false;
            }else{
                // Priority 4: LOWEST unit length (when spacer is tied)
                // Simpler repeat units (e.g., "AT" better than "ATAT")
                if(v1[0] < v2[0]){
                    return true;
                }else if(v1[0] > v2[0]){
                    return false;
                }else{
                    return false;  // Equal segmentations
                }
            }
        }
    }
}


/*
 * LTRS: Longest Tandem Repeat Search
 *
 * Counts how many consecutive times a repeat unit appears in a target sequence
 *
 * Parameters:
 *   str1_start to str1_end: The repeat unit to search for
 *   str2_start to str2_end: The target sequence to search in
 *
 * Returns: Total length of tandem repeats found (not the count, but total bases)
 *
 * Example:
 *   unit = "AT", target = "ATATGC"
 *   Returns: 4 (two copies of "AT" = 4 bases)
 */
int ltrs(std::string::iterator str1_start,std::string::iterator str1_end,std::string::iterator str2_start,std::string::iterator str2_end){

    int increment_unit = std::distance(str1_start,str1_end);  // Length of repeat unit
    int coverage = 0;  // Total bases covered by tandem repeats
    bool equal_unit=true; // assume all substrings are identical

    // Iterate through target sequence, checking if repeat unit matches
    for(auto i = str2_start; i!=str2_end;){
        for(auto j = str1_start; j!=str1_end; j++,i++){
            equal_unit = equal_unit & (*i == *j);

            if(!equal_unit){
                break;  // Unit doesn't match, stop checking
            }
        }

        if(!equal_unit){ // if is_equal_unit is false, break out.
            break;  // No more tandem repeats found
        }
        coverage += increment_unit;  // Add this repeat unit to coverage
    }
    return coverage;
}


/*
 * segmentSingle: Find optimal repeat unit segmentation for an indel
 *
 * This function tries every possible repeat unit size (1 to string length)
 * and scores each segmentation based on:
 * - How well the unit repeats internally
 * - How well it repeats in the flanking context
 * - How little "spacer" (leftover bases) remains
 *
 * Parameters:
 *   string: The indel sequence to segment
 *   context: The flanking sequence (3' context)
 *
 * Returns: Best segmentation as [unit_length, internal_reps, spacer_length, prime3_reps]
 */
std::vector<int> segmentSingle(std::string &string,std::string &context){
    // Get the total length of the indel sequence
    int string_size= string.size();

    // Create scoring matrix: one row per possible unit size, 4 columns per row
    // Column 0: unit length
    // Column 1: internal repeats (bases)
    // Column 2: spacer length (leftover bases)
    // Column 3: 3' context repeats (bases)
    std::vector<std::vector<int> > scores(string_size,std::vector<int>(4));

    int i = 0;

    // Try every possible repeat unit size (from 1 to string length)
    for(auto unit_iter = string.begin()+1; unit_iter<=string.end();++unit_iter,++i){

        // 1. Set unit length (1, 2, 3, ... string_size)
        scores[i][0]=i+1;

        // 2. Count how many times this unit repeats INSIDE the indel
        //    Uses LTRS to find tandem repeats after the first occurrence of the unit
        scores[i][1]=ltrs(string.begin(),unit_iter,unit_iter,string.end());

        // 3. Calculate spacer length (leftover bases that don't fit the repeat pattern)
        //    spacer = total_length - unit_length - internal_repeat_length
        scores[i][2]=string_size - scores[i][0] - scores[i][1];

        // 4. Count how many times this unit repeats in the 3' FLANKING CONTEXT
        //    This is crucial: indels often occur where repeats continue into flanking region
        scores[i][3]=ltrs(string.begin(),unit_iter,context.begin(),context.end());
    }

    // Sort all segmentations using priority comparison function
    // Best segmentation will be first (index 0)
    std::sort(scores.begin(),scores.end(),compareSegmentations);

    // Return the best segmentation
    return scores[0];
}
