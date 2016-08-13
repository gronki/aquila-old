Name:           aquila
Version:        607.01
Release:        1%{?dist}
Summary:        Aquila Utilities

License:        GPL 2
URL:            http://gronki.pl/aquila.tar.gz
Source0:        aquila.tar.gz

Requires:       python pyfits python-pillow dcraw numpy
BuildArch:      noarch

%description
A set of utilities useful for reduction of amateur astrophotography.

%prep
%autosetup -n aquila

%build
make

%install
rm -rf $RPM_BUILD_ROOT
%make_install prefix="/usr" datadir=%{_datadir} bindir="%{_bindir}"  libdir="%{_libdir}"  includedir="%{_includedir}"


%files
%{_bindir}/*
%{_datadir}/aquila/*


%changelog
* Fri Jul  1 2016 gronki <gronki@gmail.com>
-
