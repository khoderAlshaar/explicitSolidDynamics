import numpy as np

# Function to compute Alpha
def compute_alpha(x, y):
    return (np.degrees(np.arctan2(y, x)) + 90 + 360) % 360

# Read input data
input_file = "postProcessing/sample/0.01196965/circleIntersection_sigmaEq.xy"  # Update with the actual file name
output_file = "alphaVsSigmaEq"

data = []
with open(input_file, "r") as file:
    for line in file:
        parts = line.strip().split()
        if len(parts) == 4:
            x, y, _, sigmaEq = map(float, parts)
            alpha = compute_alpha(x, y)
            data.append((alpha, sigmaEq))

# Sort data by Alpha
data.sort()

# Write output data
with open(output_file, "w") as file:
    for alpha, sigmaEq in data:
        file.write(f"{alpha:.6f}\t{sigmaEq:.1f}\n")

print(f"Processed data written to {output_file}")
