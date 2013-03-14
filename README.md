Issues and bugs up to 2013-02-27
1. if intervals are not merged and/or intervals overlap because of the size of threshold t, than reads are added to one interval only, with the consequence, that coverage is missing on the other. This mimics underrepresented and missed areas by coverage.
WORKAROUND: merge intervals before CheKov and avoid setting threshold.
SOLUTION:  Intervals should be merged before building the TreeSet including the threshold.

