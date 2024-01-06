# Update the pnfam version number to match the git commit number

sha1=`git rev-parse --short HEAD`

date=`git log -1 --format="%ci" | cut -d" " -f1`

diff=$([ "`git diff --shortstat`" != "" ] && echo " + updates" || echo "")

cat <<HERE > version.inc
! Version
character(len=*), parameter :: version = "pnfam version ${sha1} [$date]$diff"
HERE
