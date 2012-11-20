import buildaux

# Initialize the environment with path variables, CFLAGS, and so on
bioinfo_path = '#libs/bioinfo-libs'
commons_path = '#libs/common-libs'
math_path = '#libs/math'

env = Environment(tools = ['default', 'packaging'],
                  CFLAGS = '-std=c99 -D_XOPEN_SOURCE=600 -D_GNU_SOURCE -fopenmp',
                  CPPPATH = ['#', '#src', '#include', '/usr/local/include', '/usr/include/libxml2', bioinfo_path, commons_path, math_path ],
                  LIBPATH = ['/usr/lib', '/usr/local/lib', '#libs', '#libs/common-libs/', commons_path ],
                  LIBS = ['argtable2', 'common', 'config', 'cprops', 'curl', 'gsl', 'gslcblas', 'm', 'xml2', 'z'],
                  LINKFLAGS = ['-fopenmp'])
                  
if int(ARGUMENTS.get('debug', '0')) == 1:
    debug = 1
    env['CFLAGS'] += ' -O0 -g'
else:
    debug = 0
    env['CFLAGS'] += ' -O3'

env['objects'] = []


##### Targets

# Compile dependencies
SConscript(['%s/bioformats/SConscript' % bioinfo_path,
            '%s/SConscript' % commons_path,
            '%s/SConscript' % math_path
            ], exports = ['env', 'debug'])

# Create binaries and copy them to 'bin' folder
progs = SConscript(['src/effect/SConscript',
            'src/gwas/SConscript',
            'src/vcf-tools/SConscript'
            ], exports = ['env', 'debug', 'commons_path', 'bioinfo_path', 'math_path'])

env.Install('#bin', ['hpg-variant.conf', 'vcf-info-fields.conf'])

# Create tarball
# For the packaging manager: Don't forget to point the XXX_INCLUDE_PATH and XXX_LIBRARY_PATH 
# variables to the application libraries folder!!
env.Package(NAME           = 'hpg-variant',
            VERSION        = '0.2',
            PACKAGEVERSION = 0,
            PACKAGETYPE    = 'src_targz',
            source         = env.FindSourceFiles() + env.FindHeaderFiles(progs) + 
                             [ '#libs/libargtable2.a', '#libs/libcprops.a', 
                               '#buildaux.py', '#buildvars.py', '#libs/common-libs/buildvars.py', '#libs/bioinfo-libs/buildvars.py',
                               '#COPYING', '#INSTALL' ]
            )

# Create Debian package
if 'debian' in COMMAND_LINE_TARGETS:
    SConscript("deb/SConscript", exports = ['env'] )
