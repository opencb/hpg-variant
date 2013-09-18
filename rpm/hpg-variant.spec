%define name hpg-variant
%define version 0.99.4
Name:           %{name}
Version:        %{version}
Release:        1%{?dist}
Summary:        Bioinformatics tool suite for analyzing genomic variations

License:        GPLv2
URL:            http://bioinfo.cipf.es
Source0:        %{name}-%{version}.tar.gz

BuildRequires:  gcc
BuildRequires:  glibc-devel
BuildRequires:  scons
BuildRequires:  libcurl-devel
BuildRequires:  gsl-devel
BuildRequires:  libxml2-devel
BuildRequires:  ncurses-devel
BuildRequires:  zlib-devel

Requires:       libcurl
Requires:       gsl
Requires:       libxml2
Requires:       zlib

%description
HPG Variant is a suite for the analysis of genomic data extracted from Next 
Generation Sequencing technologies. It uses parallel computing technologies 
such as CUDA and OpenMP in order to reduce the processing times.
It contains three binaries:
* hpg-var-effect retrieves the effect of genome variations.
* hpg-var-gwas conducts genomic-wide association analysis such as chi-square 
  and Fisher's exact test, and family-based analysis such as transmission 
  disequilibrium test (TDT).
* hpg-var-vcf allows one to preprocess files containing genome variations in 
  Variant Call Format.

%prep
%setup -q


%build
scons


%install
rm -rf $RPM_BUILD_ROOT

mkdir -p %{buildroot}%{_bindir}/
cp bin/hpg-var-effect %{buildroot}%{_bindir}/
cp bin/hpg-var-gwas   %{buildroot}%{_bindir}/
cp bin/hpg-var-vcf    %{buildroot}%{_bindir}/

mkdir -p %{buildroot}%{_sysconfdir}/%{name}/
cp bin/hpg-variant.conf %{buildroot}%{_sysconfdir}/%{name}/
cp bin/vcf-info-fields.conf %{buildroot}%{_sysconfdir}/%{name}/

mkdir -p %{buildroot}%{_sysconfdir}/bash_completion.d/
cp bin/bash_completion/hpg-var-effect %{buildroot}%{_sysconfdir}/bash_completion.d/
cp bin/bash_completion/hpg-var-gwas %{buildroot}%{_sysconfdir}/bash_completion.d/
cp bin/bash_completion/hpg-var-vcf %{buildroot}%{_sysconfdir}/bash_completion.d/


%files
%doc
%{_bindir}/*
%{_sysconfdir}/*


%changelog
* Wed Sep 18 2013 Cristina Yenyxe Gonzalez <cgonzalez@cipf.es> - 0.99.4
- Great reduction of memory usage in the VCF merging tool allows to merge 
  hundreds of VCF files at twice the previous speed
- VCF statistics tool retrieves statistics grouped by column (for instance, 
  population)
- Full support for Variant Call Format v4.1, including structural variants 
  and novel adjancencies with breakends
- Configuration of the Effect application simplified (entries-per-thread 
  argument removed, calculated automatically)
- Command-line autocompletion
- Arguments --version and --log-level added
- SNP ID shown in GWAS output report

* Mon Jun 10 2013 Cristina Yenyxe Gonzalez <cgonzalez@cipf.es> - 0.99.2
- VCF merging tool more tolerant to different reference alleles
- Variant effect checked for all alternate alleles of a single variant
- Memory leaks in hpg-var-effect supressed

* Wed May 07 2013 Cristina Yenyxe Gonzalez <cgonzalez@cipf.es> - 0.99.1
- New VCF filters by gene, region+type, being or not an indel, and inheritance 
  pattern (dominant/recessive)
- Multi-threaded implementation of the VCF filtering tool
- GFF/BED can be used as input for region filtering
- New VCF statistics about mendelian errors per sample, being or not an indel, 
  inheritance pattern per variant
- Possibility of saving VCF statistics to SQLite DB file
- Great reduction of memory usage, performance improvement upto 3x in VCF 
  merging tool
- VCF splitting by coverage intervals
- Miscellaneous: Some library dependencies packaged inside the application

* Mon Mar 04 2013 Cristina Yenyxe Gonzalez <cgonzalez@cipf.es> - 0.4
- Filter output uses default names, 'your_vcf_file.vcf.filtered' and
  'your_vcf_file.vcf.rejected', 'out' argument reserved for tool output
- Merge tool notifies when a sample appears more than once
- GWAS tools notify when a sample appears more than once
