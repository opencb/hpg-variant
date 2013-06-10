import SCons

def FindHeaderFiles(self, target):
    def _find_sources(tgt, src, all):
        for item in tgt:
            if SCons.Util.is_List(item):
                _find_sources(item, src, all)
            else:
                if item.abspath in all:
                    continue
                    
                all[item.abspath] = True

                if item.path.endswith('.h') and not item.path.startswith('/usr'):
                    if not item.exists():
                        item = item.srcnode()
                    src.append(item.abspath)
                else:
                    lst = item.children(scan=1)
                    _find_sources(lst, src, all)

    sources = []
    _find_sources(target, sources, {})
    return sources


def FindAllSourceFiles(self, target):
    def _find_sources(tgt, src, all):
        for item in tgt:
            if SCons.Util.is_List(item):
                _find_sources(item, src, all)
            else:
                if item.abspath in all:
                    continue
                    
                all[item.abspath] = True

                if item.path.endswith('.c'):
                    if not item.exists():
                        item = item.srcnode()
                    src.append(item.abspath)
                elif item.path.endswith('.h') and not item.path.startswith('/usr'):
                    if not item.exists():
                        item = item.srcnode()
                    src.append(item.abspath)
                else:
                    lst = item.children(scan=1)
                    _find_sources(lst, src, all)

    sources = []
    _find_sources(target, sources, {})
    return sources

        
from SCons.Script.SConscript import SConsEnvironment # just do this once
SConsEnvironment.FindHeaderFiles = FindHeaderFiles
SConsEnvironment.FindAllSourceFiles = FindAllSourceFiles
