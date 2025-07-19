# Installation

MRI2FE can be installed on macOS, Linux, and Windows systems. The package has several dependencies including ANTs, CGAL, and various Python packages.

## Prerequisites

### System Requirements

- Python 3.9, 3.10, or 3.11
- CMake 3.15 or higher
- C++ compiler (GCC, Clang, or MSVC)
- Git

### Required Software

- **ANTs (Advanced Normalization Tools)**: For image processing and registration
- **CGAL (Computational Geometry Algorithms Library)**: For mesh generation
- **vcpkg**: For managing C++ dependencies

## Installation Methods

### Option 1: Automated Installation (Recommended)

#### macOS

```bash
# Clone the repository
git clone https://github.com/turnerjennings/MRI2FE.git
cd MRI2FE

# Run the automated installation script
./install_mac.sh
```

#### Linux

```bash
# Clone the repository
git clone https://github.com/turnerjennings/MRI2FE.git
cd MRI2FE

# Run the automated installation script
./install_mac_linux.sh
```

#### Windows

```cmd
# Clone the repository
git clone https://github.com/turnerjennings/MRI2FE.git
cd MRI2FE

# Run the automated installation script
install_windows.bat
```

### Option 2: Manual Installation

1. **Install vcpkg and dependencies:**

   ```bash
   # Clone vcpkg
   git clone https://github.com/microsoft/vcpkg.git
   cd vcpkg
   ./bootstrap-vcpkg.sh  # On Windows: bootstrap-vcpkg.bat

   # Install CGAL
   ./vcpkg install cgal
   ```

2. **Install ANTs:**

   ```bash
   # Using conda (recommended)
   conda install -c conda-forge ants

   # Or build from source
   git clone https://github.com/ANTsX/ANTs.git
   cd ANTs
   mkdir build && cd build
   cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local
   make -j$(nproc)
   sudo make install
   ```

3. **Install Python dependencies:**

   ```bash
   pip install -r requirements.txt
   ```

4. **Build MRI2FE:**
   ```bash
   pip install -e .
   ```

## Verification

After installation, verify that everything is working:

```python
import MRI2FE
print("MRI2FE installed successfully!")

# Test basic functionality
from MRI2FE.utilities import COM_align
print("Core modules imported successfully!")
```

## Troubleshooting

### Common Issues

1. **CMake errors during build:**

   - Ensure CMake 3.15+ is installed
   - Check that vcpkg is properly configured
   - Verify C++ compiler is available

2. **ANTs not found:**

   - Ensure ANTs is in your PATH
   - Try installing via conda: `conda install -c conda-forge ants`

3. **CGAL dependency issues:**

   - Verify vcpkg installation
   - Check that CGAL is properly installed: `./vcpkg list | grep cgal`

4. **Python import errors:**
   - Ensure all requirements are installed: `pip install -r requirements.txt`
   - Check Python version compatibility

### Getting Help

If you encounter issues during installation:

1. Check the [GitHub Issues](https://github.com/turnerjennings/MRI2FE/issues) page
2. Review the installation logs for specific error messages
3. Ensure your system meets all prerequisites
4. Try the automated installation scripts first

## Development Installation

For developers who want to contribute to MRI2FE:

```bash
# Clone with submodules
git clone --recursive https://github.com/turnerjennings/MRI2FE.git
cd MRI2FE

# Install in development mode
pip install -e .

# Install development dependencies
pip install nox
```

## Next Steps

Once installation is complete, proceed to the [Quick Start Guide](quick-start.md) to begin using MRI2FE.
