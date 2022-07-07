FROM nfcore/base
LABEL authors="Veit Schwaemmle" \
      description="Docker image containing all requirements for tpp pipeline"

# install TPP 5.2
# copied from TPP docker file
RUN apt -y update && apt -y upgrade && apt -y install gnuplot unzip libexpat1 libexpat1-dev apache2 build-essential xsltproc nano vim && apt clean all 

RUN export PERL_MM_USE_DEFAULT=1 && cpan install CGI XML::Parser FindBin::libs JSON 

WORKDIR /usr/local/
RUN wget "http://www.tppms.org/sw/TPP5.2/TPP5.2_ubuntu.tgz" && tar xvfz TPP5.2_ubuntu.tgz && chmod -R a+rwx tpp && rm -f TPP5.2_ubuntu.tgz

WORKDIR /usr/local/tpp/conf/
RUN a2enmod rewrite && a2enmod cgid && cat httpd-2.4-tpp.conf | sed 's/#Require all granted/Require all granted/' > httpd-2.4-tpp-1.conf && perl -pi -e 's/_TPP_PORT_/10401/g' httpd-2.4-tpp-1.conf && cp -p /usr/local/tpp/conf/httpd-2.4-tpp-1.conf /etc/apache2/conf-available/httpd-tpp.conf && ln -s /etc/apache2/conf-available/httpd-tpp.conf /etc/apache2/conf-enabled/httpd-tpp.conf && perl -pi -e 's/www-data/biodocker/g' /etc/apache2/envvars

# set up conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a && rm /environment.yml

# set environment variables
ENV PERL5LIB="/usr/local/tpp/lib/perl"
ENV PATH="/usr/local/tpp/bin/:/usr/local/tpp/cgi-bin/:${PATH}"
ENV PATH="/opt/conda/envs/wombatp-tpp/bin:${PATH}"

RUN mkdir /data
WORKDIR /data