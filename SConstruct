import buildaux

# Initialize the environment with path variables, CFLAGS, and so on
bioinfo_path = '#libs/bioinfo-libs'
commons_path = '#libs/common-libs'
math_path = '#libs/math'

vars = Variables('buildvars.py')
vars.Add(PathVariable('ARGTABLE_INCLUDE_PATH', 'Path to the headers of argtable2 library', '', PathVariable.PathAccept))
vars.Add(PathVariable('ARGTABLE_LIBRARY_PATH', 'Path to the compiled argtable2 library', '', PathVariable.PathAccept))
vars.Add(PathVariable('CPROPS_INCLUDE_PATH', 'Path to the headers of cprops library', '', PathVariable.PathAccept))
vars.Add(PathVariable('CPROPS_LIBRARY_PATH', 'Path to the compiled cprops library', '', PathVariable.PathAccept))

env = Environment(variables = vars,
                  CPPPATH = ['#', '#src', '#include', '/usr/local/include', '/usr/include/libxml2', bioinfo_path, commons_path, math_path, '$ARGTABLE_INCLUDE_PATH', '$CPROPS_INCLUDE_PATH' ],
                  LIBPATH = ['/usr/lib', '/usr/local/lib', '#libs', commons_path, '$ARGTABLE_LIBRARY_PATH', '$CPROPS_LIBRARY_PATH' ],
                  LIBS = ['argtable2', 'common', 'config', 'cprops', 'curl', 'gsl', 'gslcblas', 'm', 'xml2'],
                  LINKFLAGS = ['-fopenmp'])

if int(ARGUMENTS.get('debug', '0')) == 1:
    debug = 1
    env['CFLAGS'] = '-std=c99 -D_XOPEN_SOURCE=600 -D_GNU_SOURCE -fopenmp -O0 -g'
else:
    debug = 0
    env['CFLAGS'] = '-std=c99 -D_XOPEN_SOURCE=600 -D_GNU_SOURCE -fopenmp -O3'

env['objects'] = []

# Targets

SConscript(['%s/bioformats/SConscript' % bioinfo_path,
            '%s/SConscript' % commons_path,
            '%s/SConscript' % math_path
            ], exports = ['env', 'debug'])

SConscript(['src/effect/SConscript',
            'src/gwas/SConscript',
            'src/vcf-tools/SConscript'
            ], exports = ['env', 'debug', 'commons_path', 'bioinfo_path', 'math_path'])

env.Install('#bin', ['hpg-variant.cfg', 'vcf-info-fields.cfg'])
