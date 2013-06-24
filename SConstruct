import buildaux

# Initialize the environment with path variables, CFLAGS, and so on
bioinfo_path = '#libs/bioinfo-libs'
commons_path = '#libs/common-libs'
math_path = '#libs/math'

env = Environment(tools = ['default', 'packaging'],
                  CFLAGS = '-std=c99 -D_XOPEN_SOURCE=600 -D_GNU_SOURCE -msse4.2 -fopenmp',
                  CPPPATH = ['#', '#src', '#include', bioinfo_path, commons_path, math_path, 
                             '/usr/include', '/usr/local/include', '/usr/include/libxml2', '/usr/lib/openmpi/include'],
                  LIBPATH = [commons_path, bioinfo_path, '/usr/lib', '/usr/local/lib'],
                  LIBS = ['common', 'bioinfo', 'curl', 'dl', 'gsl', 'gslcblas', 'm', 'mpi', 'xml2', 'z'],
                  LINKFLAGS = ['-fopenmp'])
                  
if int(ARGUMENTS.get('debug', '0')) == 1:
    debug = 1
    env['CFLAGS'] += ' -O0 -g'
else:
    debug = 0
    env['CFLAGS'] += ' -O3 -g'

env['objects'] = []

# bioinfo-libs compilation
formats = [ 'family', 'features', 'gff', 'ped', 'vcf' ]
aligners = []
compiler = 'gcc'

##### Targets

# Compile dependencies
SConscript(['%s/SConscript' % bioinfo_path,
            '%s/SConscript' % commons_path,
            '%s/SConscript' % math_path
            ], exports = ['env', 'debug', 'formats', 'aligners', 'compiler'])

# Create binaries and copy them to 'bin' folder
progs = SConscript(['src/effect/SConscript',
            'src/gwas/SConscript',
            'src/vcf-tools/SConscript'
            ], exports = ['env', 'debug', 'commons_path', 'bioinfo_path', 'math_path'])

inst = env.Install('#bin', ['hpg-variant.conf', 'vcf-info-fields.conf'] + list(progs))

# Run tests
t = SConscript("test/SConscript", exports = ['env', 'debug', 'commons_path', 'bioinfo_path', 'math_path'] )

# Create tarball
# For the packaging manager: Don't forget to point the XXX_INCLUDE_PATH and XXX_LIBRARY_PATH 
# variables to the application libraries folder!!
tb = env.Package(NAME          = 'hpg-variant',
                VERSION        = '0.99.2',
                PACKAGEVERSION = 0,
                PACKAGETYPE    = 'src_targz',
                source         = env.FindSourceFiles() + env.FindHeaderFiles(progs) + 
                             [ '#buildaux.py', '#libs/bioinfo-libs/buildvars.py',
                               '#deb/SConscript', '#rpm/SConscript', '#rpm/hpg-variant.spec',
                               '#COPYING', '#INSTALL' ] )
Alias('tarball', tb)
  
# Create Debian package
deb = SConscript("deb/SConscript", exports = ['env'] )
Alias('debian', deb)

#Create Fedora package
# TODO Should not be necessary to set a condition here! And 'tarball' must be run by hand :(
if 'fedora' in COMMAND_LINE_TARGETS:
    fed = SConscript("rpm/SConscript", exports = ['env'] )
    Depends(fed, tb)
    Alias('fedora', fed)


# By default, create only the executables and install them
Default(progs, inst)
