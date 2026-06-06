# Agent Guidelines for FitHBonds

## Terminal Command Behavior

- **Never use `exit` in run_command** - causes "Exit code could not be determined"
- Good: `echo "test"` or `ls *.py | head -5`
- Bad: `echo "test"; exit 0` or `ls *.py; exit 0`
