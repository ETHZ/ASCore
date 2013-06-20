REMOTE=git://github.com/ETHZ
PACKNAME=$1
VERSION=$2
PACKNAME_GIT=$(echo ${PACKNAME} | sed 's,/,-,').git

echo Checking out from the ETHZ github area the package ${PACKNAME} version ${VERSION}
git clone ${REMOTE}/${PACKNAME_GIT} ${PACKNAME}
cd ${PACKNAME}
git checkout ${VERSION}
cd ${OLDPWD}
