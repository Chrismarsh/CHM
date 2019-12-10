#/bin/bash

# This runs the conan build missing command when running under travis-ci.
# The build has too large of a log output, so we need to redirect to dev null
# however if travis-ci hasn't seen stdout for   10min, it'll terminate the build.
# This busy loops waiting for the build to finish while periodically updating stdout

(conan install . --build missing) &> /dev/null &
PID=$!
while [ -d /proc/$PID ]
do
    echo "Building..."
    sleep 5
done