"""Python helper tools for CMS FWLite

benedikt.hegner@cern.ch

"""
import re
import ROOT
import exceptions
### define tab completion
try:
  import readline, cmscompleter
  readline.parse_and_bind('tab: complete')
except:
  print 'WARNING: Could not load tab completion'


# for adding iterators at runtime
import iterators


### workaround iterator generators for ROOT classes
def all(container):

  # loop over ROOT::TTree and similar
  if hasattr(container,'GetEntries'):
    try:
      entries = container.GetEntries()
      for entry in xrange(entries):
        yield entry
    except:
        raise cmserror("Looping of %s failed" %container) 

  # loop over std::vectors and similar
  elif hasattr(container, 'size'):
    try:
      entries = container.size()
      for entry in xrange(entries):
        yield container[entry]
    except:
      pass

 # loop over containers with begin and end iterators
def loop(begin, end):
    """Convert a pair of C++ iterators into a python generator"""
    while (begin != end):
        yield begin.__deref__()  #*b
        begin.__preinc__()       #++b 

### auto branch types (Chris Jones)
def createBranchBuffer(branch):
    reColons = re.compile(r'::')
    reCloseTemplate =re.compile(r'>')
    reOpenTemplate =re.compile(r'<')
    branchType = ROOT.branchToClass(branch)
    buffer = eval ('ROOT.'+reColons.sub(".",reOpenTemplate.sub("(ROOT.",reCloseTemplate.sub(")",branchType.GetName())))+'()')
    if( branch.GetName()[-1] != '.'):
        branch.SetAddress(buffer)
    else:
        branch.SetAddress(ROOT.AddressOf(buffer))
    return buffer


class EventTree(object):
      def __init__(self,obj):
          if isinstance(obj, ROOT.TTree):
              self._tree = obj
          elif isinstance(obj, ROOT.TFile):
              self._tree = obj.Get("Events")
          elif isinstance(obj, str):
              self._tree = ROOT.TFile.Open(obj).Get("Events")
          else:
              raise cmserror("EventTree accepts only TTrees, TFiles and filenames")
          self._usedBranches = dict()
          self._index = -1
          self._aliases = self._tree.GetListOfAliases()
      def branch(self,name):
          # support for aliases
          alias = self._tree.GetAlias(name)
          if alias != '': name = alias 
          # access the branch in ttree
          if name in self._usedBranches:
              return self._usedBranches[name]
          self._usedBranches[name]=EventBranch(self,name)
          return self._usedBranches[name]
      def cppCode(self, name):
          """C++ code for accessing the product inside the full framework"""
          alias = self._tree.GetAlias(name)
          if alias != '': name = alias
          tmpBranch = self._tree.GetBranch(name)
          typeString = ROOT.branchToClass(tmpBranch).GetName()
          nameParts = name.split("_")
          if nameParts[2] == "":
              cppCode = 'edm::Handle<%s > dummy;\nevent.getByLabel("%s", dummy);'\
                        %(typeString, nameParts[1])
          else:
              cppCode = 'edm::Handle<%s > dummy;\nevent.getByLabel("%s", "%s", dummy);'\
                        %(typeString, nameParts[1], nameParts[2])
          return cppCode
      def getListOfAliases(self):
          return self._aliases
      def index(self):
          return self._index
      def tree(self):
          return self._tree
      def __setBranchIndicies(self):
          for branch in self._usedBranches.itervalues():
              branch.setIndex(self._index)
      def __getattr__(self, name):
          return self.branch()
      def __getitem__(self,key):
          if key <0 or key > self._tree.GetEntries():
              raise IndexError
          self._index = key
          self.__setBranchIndicies()
          return self
      def __iter__(self):
          for entry in xrange(self._tree.GetEntries()):
              self._index = entry
              self.__setBranchIndicies()
              self._tree.GetEntry(self._index,0)
              yield Event(self)    # TODO: don't return a new object but update the old one 
              

class Event(object):
    def __init__(self, eventTree):
        self._eventTree = eventTree

    def getProduct(name):
        return iterators.addIterator(self._eventTree.branch(name)())

    def __getattr__(self, name):
        return iterators.addIterator(self._eventTree.branch(name)())


class EventBranch(object):
    def __init__(self,parent,name):
        self._branch = parent.tree().GetBranch(name)
        if self._branch == None:
            raise cmserror("Unknown branch "+name)
        self._buffer = createBranchBuffer(self._branch)
        self._index = parent.index()
        self._readData = False
    def setIndex(self,index):
        self._index = index
        self._readData = False
    def __readData(self):
        self._branch.GetEntry(self._index)
        self._readData = True

    # replace this by __getattr__ to allow branch.attr instead of branch().attr
    def __call__(self):
        if not self._readData:
            self.__readData()
        return self._buffer


class cmserror(exceptions.StandardError):
    def __init__(self, message):
          length = len(message)+7   #7=len("ERROR: ")
          print "="*length
          print "ERROR:", message
          print "="*length
