export BASEDIR=$PWD

export PATH=$BASEDIR/install/bin:$PATH

echo "*********************** mkedm to compile *******************"
alias mkedm="cd $BASEDIR/build ; make -j8 install ; cd -"
