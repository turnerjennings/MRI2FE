# Contributing to MRI2FE

## Getting Started

### Prerequisites

Before you begin, ensure you have the following installed:

- Python 3.8 or higher
- Git
- pip

### Setting Up the Development Environment

1. **Fork the repository** on GitHub by clicking the "Fork" button.

2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/turnerjennings/MRI2FE.git
   cd MRI2FE
   ```

3. **Create a virtual environment**:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

4. **Install the package in development mode**:
   ```bash
   pip install -e ".[dev]"
   ```

5. **Install nox to automate build and test functions**:
    ```bash
    pip install nox
    ```

## Development Workflow

### Making Changes

1. Make your changes in the appropriate files
2. Write or update tests for your changes
3. Update documentation if necessary
4. Ensure your code follows our style guidelines

### Running Tests

Before submitting your changes, run the nox suite to test, lint, and format the code:

```bash
#run all test and lint
nox

#run only tests
nox -s test
nox -s cpptest

#run lint nd format
nox -s lint
nox -s format+

```

## Types of Contributions

### Bug Reports

When reporting bugs, please include:

- A clear, descriptive title
- Steps to reproduce the issue
- Expected vs actual behavior
- Your environment details (Python version, OS, package version)
- Minimal code example that demonstrates the issue

### Feature Requests

For new features:

- Explain the motivation and use case
- Provide a detailed description of the proposed functionality
- Consider backward compatibility
- Include examples of how the feature would be used

### Code Contributions

#### Pull Request Process

1. **Update documentation** for any new features or API changes
2. **Add tests** that cover your changes
3. **Ensure all tests pass** and code quality checks succeed
4. **Create a pull request** with a clear title and description

### Documentation

Documentation improvements are always welcome! This includes:

- Fixing typos or clarifying existing documentation
- Adding examples or tutorials
- Improving API documentation
- Translating documentation

### Project Structure

```
your-package-name/
├── docs/                       # Documentation source
├── include/                    # cpp header files
├── src/
    ├── cpp/                    # cpp code files
    ├── MRI2FE/                 # python package files    
├── test/                       # Test files
├── noxfile.py                  # nox session configuration
├── pyproject.toml              # Project configuration
├── install_windows.bat         # Windows installation script
├── install_linux.sh            # Linux installation script
├── install_mac.sh              # MacOS installation script
├── pyproject.toml              # Project configuration
├── README.md                   # Project overview
└── CONTRIBUTING.md             # This file
```
