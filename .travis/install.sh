#!/bin/bash

set -e
set -x

if [[ "$(uname -s)" == 'Darwin' ]]; then
    brew update || brew update
    brew outdated pyenv || brew upgrade pyenv
    brew install pyenv-virtualenv
    brew install cmake || true

    if which pyenv > /dev/null; then
        eval "$(pyenv init -)"
    fi

    pyenv install 2.7.16
    pyenv virtualenv 2.7.16 conan
    pyenv rehash
    pyenv activate conan
fi

pip install conan --upgrade
pip install conan_package_tools

conan profile new default --detect
conan remote add bincrafters "https://api.bintray.com/conan/bincrafters/public-conan"
conan remote add CHM "https://api.bintray.com/conan/chrismarsh/CHM"
conan profile update settings.compiler.cppstd=14 default
conan profile update settings.compiler.libcxx=libstdc++11 default

conan user
