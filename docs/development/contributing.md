# Contributing to MRI2FE

Thank you for your interest in contributing to MRI2FE! This guide will help you get started with development and contributing to the project.

## Development Setup

### Prerequisites

1. **Python Environment**: Python 3.9, 3.10, or 3.11
2. **Development Tools**: Git, CMake, C++ compiler
3. **Dependencies**: ANTs, CGAL, vcpkg

### Setting Up Development Environment

1. **Clone the repository with submodules:**

   ```bash
   git clone --recursive https://github.com/turnerjennings/MRI2FE.git
   cd MRI2FE
   ```

2. **Install in development mode:**

   ```bash
   pip install -e .
   ```

3. **Install development dependencies:**

   ```bash
   pip install nox ruff black pytest
   ```

4. **Verify installation:**
   ```bash
   python -c "import MRI2FE; print('Development setup complete!')"
   ```

## Code Style and Standards

### Python Code Style

MRI2FE follows PEP 8 with some modifications:

- **Line length**: 88 characters (Black default)
- **Docstrings**: Google style
- **Type hints**: Required for public functions
- **Import organization**: Standard library, third-party, local

### Running Code Quality Checks

```bash
# Format code with Black
nox -s format

# Lint code with Ruff
nox -s lint

# Run tests
nox -s test

# Run all checks
nox -s lint-format
```

### Pre-commit Hooks

Install pre-commit hooks to automatically format and lint code:

```bash
pip install pre-commit
pre-commit install
```

## Testing

### Running Tests

```bash
# Run all tests
pytest

# Run specific test file
pytest test/test_utilities.py

# Run with coverage
pytest --cov=MRI2FE

# Run tests in parallel
pytest -n auto
```

### Writing Tests

Tests should be written using pytest. Follow these guidelines:

1. **Test file naming**: `test_<module_name>.py`
2. **Test function naming**: `test_<function_name>_<scenario>`
3. **Use fixtures** for common setup
4. **Test both success and failure cases**
5. **Include edge cases and error conditions**

Example test:

```python
import pytest
import numpy as np
from MRI2FE.utilities import COM_align

def test_COM_align_basic():
    """Test basic COM alignment functionality."""

    # Create test data
    fixed = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    moving = np.array([[1, 1, 1], [2, 1, 1], [1, 2, 1]])

    # Perform alignment
    result = COM_align(fixed, moving)

    # Check results
    assert result.shape == moving.shape
    assert not np.allclose(result, moving)  # Should be different after alignment

def test_COM_align_invalid_input():
    """Test COM alignment with invalid inputs."""

    with pytest.raises(ValueError, match="fixed cannot be None"):
        COM_align(None, np.array([[0, 0, 0]]))

    with pytest.raises(ValueError, match="moving cannot be None"):
        COM_align(np.array([[0, 0, 0]]), None)
```

### Test Data

- Store test data in `test/test_data/`
- Use small, representative datasets
- Include both valid and invalid test cases
- Document test data sources and formats

## Documentation

### Docstring Standards

All public functions and classes must have docstrings in Google style:

```python
def function_name(param1: str, param2: int = 1) -> bool:
    """Brief description of function.

    Longer description if needed.

    Args:
        param1: Description of param1
        param2: Description of param2. Defaults to 1.

    Raises:
        ValueError: When param1 is invalid
        TypeError: When param2 is wrong type

    Returns:
        Description of return value

    Example:
        >>> function_name("test", 2)
        True
    """
    pass
```

### Updating Documentation

1. **Update docstrings** in source code
2. **Update user guide** pages as needed
3. **Add examples** for new features
4. **Update API reference** (automatic with mkdocstrings)

### Building Documentation

```bash
# Build documentation locally
mkdocs serve

# Build for deployment
mkdocs build
```

## Development Workflow

### 1. Create a Feature Branch

```bash
git checkout -b feature/your-feature-name
```

### 2. Make Changes

- Write code following style guidelines
- Add tests for new functionality
- Update documentation
- Update changelog if needed

### 3. Test Your Changes

```bash
# Run all quality checks
nox -s lint-format

# Run tests
nox -s test

# Build documentation
mkdocs build
```

### 4. Commit Your Changes

```bash
git add .
git commit -m "feat: add new feature description

- Detailed description of changes
- Any breaking changes
- Related issues"
```

### 5. Push and Create Pull Request

```bash
git push origin feature/your-feature-name
```

## Pull Request Guidelines

### Before Submitting

1. **Ensure all tests pass**
2. **Update documentation**
3. **Add changelog entry**
4. **Check code coverage**
5. **Review your own code**

### Pull Request Template

```markdown
## Description

Brief description of changes

## Type of Change

- [ ] Bug fix
- [ ] New feature
- [ ] Breaking change
- [ ] Documentation update

## Testing

- [ ] Tests added/updated
- [ ] All tests pass
- [ ] Manual testing completed

## Documentation

- [ ] Docstrings updated
- [ ] User guide updated
- [ ] API reference updated

## Checklist

- [ ] Code follows style guidelines
- [ ] Self-review completed
- [ ] Changelog updated
- [ ] No breaking changes (or documented)
```

## Release Process

### Versioning

MRI2FE follows [Semantic Versioning](https://semver.org/):

- **MAJOR**: Breaking changes
- **MINOR**: New features, backward compatible
- **PATCH**: Bug fixes, backward compatible

### Creating a Release

1. **Update version** in `pyproject.toml`
2. **Update changelog**
3. **Create release branch**
4. **Run full test suite**
5. **Build and test documentation**
6. **Create GitHub release**
7. **Publish to PyPI**

## Common Development Tasks

### Adding a New Function

1. **Add function to appropriate module**
2. **Add type hints and docstring**
3. **Write tests**
4. **Update documentation**

### Adding a New Module

1. **Create module file**
2. **Add `__init__.py` imports**
3. **Write comprehensive tests**
4. **Add to documentation**

### Fixing a Bug

1. **Create test that reproduces bug**
2. **Fix the bug**
3. **Ensure test passes**
4. **Add regression test**

## Getting Help

### Development Resources

- **Issues**: Report bugs and request features
- **Discussions**: Ask questions and share ideas
- **Code Review**: Get feedback on your changes

### Communication

- **Be respectful** and inclusive
- **Provide context** when asking questions
- **Use clear, descriptive language**
- **Follow up** on discussions

## Code of Conduct

Please read and follow our [Code of Conduct](CODE_OF_CONDUCT.md). We are committed to providing a welcoming and inclusive environment for all contributors.

## License

By contributing to MRI2FE, you agree that your contributions will be licensed under the MIT License.
