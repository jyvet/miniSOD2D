#!/bin/bash

# Output file to store total execution time
output_file="execution_times.txt"

# Clear the output file if it already exists
> "$output_file"

# Start the timer
start_time=$(date +%s.%N)

# Loop 100 times
for i in {1..10}
do
    ./src/app/miniSOD
done

# End the timer
end_time=$(date +%s.%N)

# Calculate the total elapsed time
total_time=$(echo "$end_time - $start_time" | bc)

# Write the total time to the output file
echo "Total time for 100 runs: $total_time seconds" > "$output_file"

echo "Total execution time recorded in $output_file"
