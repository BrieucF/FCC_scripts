export EOS_MGM_URL=root://eospublic.cern.ch
username=mgarciam
firstletter=${username:0:1}
userid="$(id -u $username)"
usergid="$(id -g $username)"
storagepath=/eos/experiment/fcc/users/$firstletter/$username/
echo mkdir -p $storagepath
echo eos chown $userid:$usergid $storagepath 
echo eos quota set -u $username -v 2T -p $storagepath
