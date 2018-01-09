rm -rf nestml_build
mkdir nestml_build

NESTML_INSTALL_DIR="${NESTML_INSTALL_DIR:-$HOME/workspace/3rdparty/nestml/target}"
java -jar "$NESTML_INSTALL_DIR"/nestml.jar custom_hh_model/ --target nestml_build/

cd nestml_build
cmake -Dwith-nest=$NEST_INSTALL_DIR/bin/nest-config .
make all
make install

