
env = Environment()

# Check dependency libraries: cprops
conf = Configure(env)

if not conf.CheckLib('argtable2'):
    print 'argtable2 library not found!'
    Exit(1)

if not conf.CheckLib('config'):
    print 'config library not found!'
    Exit(1)

if not conf.CheckLib('cprops'):
    print 'cprops library not found!'
    Exit(1)

if not conf.CheckLib('curl'):
    print 'cURL library not found!'
    Exit(1)

if not conf.CheckLib('gsl'):
    print 'GSL library not found!'
    Exit(1)

if not conf.CheckLib('xml2'):
    print 'xml2 library not found!'
    Exit(1)
    
env = conf.Finish()


# Extra flags
bioinfo_path = '#libs/bioinfo-libs'
commons_path = '#libs/common-libs'
math_path = '#libs/math'

if int(ARGUMENTS.get('debug', '0')) == 1:
    debug = 1
    env['CFLAGS'] = '-std=c99 -D_XOPEN_SOURCE=600 -D_GNU_SOURCE -fopenmp -O0 -g'
else:
    debug = 0
    env['CFLAGS'] = '-std=c99 -D_XOPEN_SOURCE=600 -D_GNU_SOURCE -fopenmp -O3'
    
env['CPPPATH'] = ['#', '#src', '#include', bioinfo_path, commons_path, math_path]
env['LIBPATH'] = ['/usr/lib', '/usr/local/lib', '#libs', commons_path]
env['LIBS']   += ['common']
env['LINKFLAGS'] = '-fopenmp'

env.ParseConfig("pkg-config argtable2 --cflags --libs")
env.ParseConfig("pkg-config libconfig --cflags --libs")
env.ParseConfig("pkg-config libcurl --cflags --libs")
env.ParseConfig("pkg-config gsl --cflags --libs")
env.ParseConfig("pkg-config libxml-2.0 --cflags --libs")
env.ParseConfig("pkg-config libcurl --cflags --libs")

print 'CPPPATH = %s' % env['CPPPATH']
print 'LIBPATH = %s' % env['LIBPATH']
print 'LIBS = %s' % env['LIBS']


## Targets

SConscript(['%s/bioformats/SConscript' % bioinfo_path,
            '%s/SConscript' % commons_path,
            '%s/SConscript' % math_path
            ], exports = ['env', 'debug'])

SConscript(['src/effect/SConscript',
            'src/gwas/SConscript',
            'src/vcf-tools/SConscript'
            ], exports = ['env', 'debug', 'commons_path', 'bioinfo_path', 'math_path'])

