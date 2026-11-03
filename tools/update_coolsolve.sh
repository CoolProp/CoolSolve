#!/bin/bash
#
# update_coolsolve.sh -- update and restart the CoolSolve online server.
#
#   *** RUN AS ROOT ***   sudo /home/coolsolve/CoolSolve/tools/update_coolsolve.sh
#
# Root is needed to restart the services; the git pull, the GUI build and the
# compilation run as the unprivileged user that owns the checkout ($USER below).
# The layout is the one of docs/deployment_ubuntu_apache.md (checkout in
# /home/coolsolve/CoolSolve, build directory build/, systemd service coolsolve).
#
# What gets updated:
#   - CoolSolve: `git pull` of the current branch.
#   - GUI: rebuilt from the lockfile (npm ci) and embedded in the binary.
#   - CoolSolve Library: `cmake ..` re-fetches its `main` branch, so the binary
#     embeds the latest library (when GitHub is unreachable cmake warns and keeps
#     the copy already fetched).
#   - CoolProp: the commit pinned in CMakeLists.txt (COOLSOLVE_COOLPROP_PINNED_TAG).
#     `-U COOLSOLVE_COOLPROP_TAG` drops the tag cached by an earlier configure, so a
#     pin bump pulled from git is followed (a bump recompiles CoolProp: long).
#
# Check the result in the last lines it prints (`coolsolve --version`).

set -e

PROJECT_ROOT="/home/coolsolve/CoolSolve"
USER="coolsolve"
SERVICE_NAME="coolsolve.service"
APACHE_SERVICE="apache2.service"
# Few jobs: embedding the model library needs about 1 GB of RAM per compile job
JOBS=4

if [ "$(id -u)" -ne 0 ]; then
    echo "ERROR: this script must be run as root (sudo $0)" >&2
    exit 1
fi

echo "Starting CoolSolve update process..."

# 1. Pull latest changes from git
echo "Pulling latest changes from git..."
sudo -u $USER -i bash -c "cd $PROJECT_ROOT && git pull"

# 2. Rebuild the GUI
echo "Rebuilding the GUI..."
sudo -u $USER -i bash -c "cd $PROJECT_ROOT/gui && npm ci && npm run build"

# 3. Recompile the backend (only the coolsolve binary, not the test suite)
echo "Recompiling the backend..."
sudo -u $USER -i bash -c "mkdir -p $PROJECT_ROOT/build && cd $PROJECT_ROOT/build && cmake -U COOLSOLVE_COOLPROP_TAG -DCOOLSOLVE_BUILD_GUI=ON .. && make -j$JOBS coolsolve"

# 4. Restart the CoolSolve service
echo "Restarting CoolSolve service..."
systemctl restart $SERVICE_NAME

# 5. Reload Apache
echo "Reloading Apache..."
systemctl reload $APACHE_SERVICE

echo "Update complete! Now serving:"
sudo -u $USER -i bash -c "$PROJECT_ROOT/build/coolsolve --version"
