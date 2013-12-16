import buildaux

# Initialize the environment with path variables, CFLAGS, and so on
bioinfo_path = '#lib/bioinfo-libs'
commons_path = '#lib/common-libs'
math_path = '#lib/math'

compiler = ARGUMENTS.get('compiler', 'gcc')
build_tools = [ 'default', 'packaging' ]
if compiler == 'icc':
    build_tools.append('intelc')

env = Environment(tools = build_tools,
                  CFLAGS = '-std=c99 -D_XOPEN_SOURCE=600 -D_GNU_SOURCE -fopenmp -Wuninitialized -Wmissing-braces',
                  CPPPATH = ['#', '#src', '#include', bioinfo_path, commons_path, math_path, '/usr/include', '/usr/local/include', '/usr/include/libxml2'],
                  LIBPATH = [commons_path, bioinfo_path, '/usr/lib', '/usr/local/lib'],
                  LIBS = ['common', 'bioinfo', 'curl', 'dl', 'gsl', 'gslcblas', 'm', 'xml2', 'z'],
                  LINKFLAGS = ['-fopenmp'])
                  
if int(ARGUMENTS.get('debug', '0')) == 1:
    debug = 1
    env['CFLAGS'] += ' -O0 -g'
else:
    debug = 0
    env['CFLAGS'] += ' -O3 -g'

env['objects'] = []

##### Targets

# Compile dependencies
SConscript(['%s/SConscript' % bioinfo_path,
            '%s/SConscript' % commons_path,
            '%s/SConscript' % math_path
            ], exports = ['env', 'debug', 'compiler'])

# Create binaries and copy them to 'bin' folder
progs = SConscript(['src/effect/SConscript',
            'src/gwas/SConscript',
            'src/vcf-tools/SConscript'
            ], exports = ['env', 'debug', 'commons_path', 'bioinfo_path', 'math_path'])

inst1 = env.Install('#bin', [Glob('etc/hpg-variant/*.conf')] + list(progs))
inst2 = env.Install('#bin/bash_completion', [Glob('etc/bash_completion.d/*')])


# Run tests
t = SConscript("test/SConscript", exports = ['env', 'debug', 'commons_path', 'bioinfo_path', 'math_path'] )

# Create tarball
# For the packaging manager: Don't forget to point the XXX_INCLUDE_PATH and XXX_LIBRARY_PATH 
# variables to the application libraries folder!!
tb = env.Package(NAME          = 'hpg-variant',
                VERSION        = '1.0',
                PACKAGEVERSION = 0,
                PACKAGETYPE    = 'src_targz',
                source         = env.FindSourceFiles() + env.FindHeaderFiles(progs) + 
                             [ '#buildaux.py', '#lib/bioinfo-libs/buildvars.py',
                               '#deb/SConscript', '#rpm/SConscript', '#rpm/hpg-variant.spec',
                               Glob('etc/hpg-variant/*.conf'), Glob('etc/bash_completion.d/*'),
                               '#COPYING', '#README' ] )
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
Default(progs, inst1, inst2)
