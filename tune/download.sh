wget -qO- http://www.fftw.org/fftw-3.3.9.tar.gz | tar -xvz -C ./dependencies/ 
cd dependencies/fftw-3.3.9 && ./configure --prefix  $(pwd)../../build && make && make install

git clone https://github.com/sciplot/sciplot --recursive dependencies/sciplot
cd dependencies/sciplot  && cmake . -DCMAKE_INSTALL_PREFIX=$(pwd)../../build && cmake --build . --target install

