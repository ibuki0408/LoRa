# Use the official Julia image as the base
FROM julia:1.10-bookworm

# Set the working directory in the container
WORKDIR /app

# Copy the Project.toml to the container
# (Manifest.toml is optional, but Project.toml is required)
COPY Project.toml ./

# Install dependencies
# We use Pkg.instantiate() to install the exact versions if possible,
# or Pkg.resolve() if Manifest.toml was missing.
RUN julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# Copy the rest of the application code
# (Note: In development/server-run, we often overwrite this with a volume mount)
COPY . .

# Set the default command to run the simulation script
# The script itself handles parallelization based on available cores.
CMD ["julia", "--project", "scripts/run_prop_per.jl"]
