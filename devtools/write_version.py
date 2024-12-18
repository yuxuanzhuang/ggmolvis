import ggmolvis

# Extract the version
version = ggmolvis.__version__

# Write to a file
with open("VERSION.txt", "w") as version_file:
    version_file.write(version + "\n")