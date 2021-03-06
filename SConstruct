import os
import buildaux

# Initialize the environment with path variables, CFLAGS, and so on
hpglib_path = '#lib/c/'
third_party_path = '#/lib/third_party'
#third_party_cram_path = '#/lib/third_party/htslib/cram'
third_party_hts_path = '#/lib/third_party/htslib'
third_party_samtools_path = '#/lib/third_party/samtools'

#bioinfo_path = '#lib/bioinfo-libs'
#commons_path = '#lib/common-libs'
#math_path = '#lib/math'

compiler = ARGUMENTS.get('compiler', 'gnu')
build_tools = [ 'default', 'packaging' ]
if compiler == 'intel':
    build_tools.append('intelc')

env = Environment(tools = build_tools,
                  CFLAGS = '-std=c99 -D_XOPEN_SOURCE=600 -D_GNU_SOURCE -msse4.2 -fopenmp -Wuninitialized -Wmissing-braces',
                  CPPPATH = ['#', '#src', '#include', hpglib_path + 'src', third_party_path, third_party_samtools_path, third_party_hts_path,
                             '/usr/include', '/usr/local/include', '/usr/include/libxml2', '/usr/lib/openmpi/include'],
                  LIBPATH = [hpglib_path + 'build', third_party_hts_path, third_party_samtools_path, '/usr/lib', '/usr/local/lib'],
                  LIBS = ['curl', 'dl', 'gsl', 'gslcblas', 'm', 'xml2', 'z'],
                  LINKFLAGS = ['-fopenmp'])


if os.environ.has_key('CPATH'):
    for dir in os.getenv('CPATH').split(':'):
        env.Append(CPPPATH=[dir])

if os.environ.has_key('LIBRARY_PATH'):
    for dir in os.getenv('LIBRARY_PATH').split(':'):
        env.Append(LIBPATH=[dir])


mode = 'single'
if ARGUMENTS.get('mode', 'single') == 'mpi':
    env['CFLAGS'] += ' -D_USE_MPI'
    env['LIBS'] += ['mpi']
    mode = 'mpi'


if int(ARGUMENTS.get('debug', '0')) == 1:
    debug = 1
    env['CFLAGS'] += ' -O0 -g'
else:
    debug = 0
    env['CFLAGS'] += ' -O3 -g'

env['objects'] = []

##### Targets

# Compile dependencies
SConscript(['lib/SConstruct'])
#SConscript(['%s/SConscript' % bioinfo_path,
#'%s/SConscript' % commons_path,
#'%s/SConscript' % math_path
#], exports = ['env', 'debug', 'compiler'])

# Create binaries and copy them to 'bin' folder
progs = SConscript(['src/effect/SConscript',
            'src/gwas/SConscript',
            'src/vcf-tools/SConscript'
            ], exports = ['env', 'debug', 'mode', 'hpglib_path', 'third_party_hts_path', 'third_party_samtools_path'])

inst1 = env.Install('#bin', [Glob('etc/hpg-variant/*.conf')] + list(progs))
inst2 = env.Install('#bin/bash_completion', [Glob('etc/bash_completion.d/*')])


# Run tests
t = SConscript("test/SConscript", exports = ['env', 'debug', 'hpglib_path', 'third_party_hts_path', 'third_party_samtools_path'] )

# Create tarball
# For the packaging manager: Don't forget to point the XXX_INCLUDE_PATH and XXX_LIBRARY_PATH 
# variables to the application libraries folder!!
tb = env.Package(NAME          = 'hpg-variant',
                VERSION        = '2.0',
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

