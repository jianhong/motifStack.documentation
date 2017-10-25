This is a helper documentation of motifStack.

[Docker](https://docs.docker.com/) container allows software to be packaged into containers which can be run in any platform. Users can download the motifStack docker image using the following code snippet and render the rmarkdown files.

In Docker jianhong/motifstack, Matalign and phylip/neighbor are installed at /usr/bin. In Docker jianhong/motifstack_1.21.6, git, texlive and pandoc are also installed.

## how to run
<pre>
## in windows, open docker terminal
cd ~ ## in windows, please try cd c:\\ Users\\ username
docker pull jianhong/motifstack:latest
mkdir tmp4motifstack ## this will be the share folder for your host and container.
docker run -ti --rm -v ${PWD}/tmp4motifstack:/volume/data jianhong/motifstack:latest bash
  In motifstack:latest docker
    1  cd /volume/data
    2  git clone https://github.com/jianhong/motifStack.documentation.git
    3  cd motifStack.documentation/
    4  cp /usr/bin/matalign app/matalign-v4a
    5  cp /usr/bin/phylip/neighbor app/neighbor.app/Contents/MacOS/neighbor
    6  R cmd -e "rmarkdown::render('suppFigure2.Rmd')"
    7  R cmd -e "rmarkdown::render('suppFigure6.Rmd')"
</pre>