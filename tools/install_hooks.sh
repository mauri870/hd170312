#!/bin/sh
# Install the hooks that run tools/scan_names.py before each commit and on each commit message.
set -e
root=$(git rev-parse --show-toplevel)
printf '#!/bin/sh\nexec python3 "%s/tools/scan_names.py" --staged\n' "$root" > "$root/.git/hooks/pre-commit"
printf '#!/bin/sh\nexec python3 "%s/tools/scan_names.py" --message "$1"\n' "$root" > "$root/.git/hooks/commit-msg"
chmod +x "$root/.git/hooks/pre-commit" "$root/.git/hooks/commit-msg"
