package confdb.data;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.Collections;


/**
 * Configuration
 * -------------
 * @author Philipp Schieferdecker
 *
 * Description of a CMSSW job configuration.
 */
public class Configuration implements IConfiguration
{
    //
    // member data
    //

    /** configuration information */
    private ConfigInfo configInfo = null;
    
    /** current software release */
    private SoftwareRelease release = null;
    
    /** has the configuration changed since the last 'save' operation? */
    private boolean hasChanged = false;
    
    /** list of globale parameter sets */
    private ArrayList<PSetParameter>    psets = null;

    /** list of EDSources */
    private ArrayList<EDSourceInstance> edsources = null;

    /** list of ESSources */
    private ArrayList<ESSourceInstance> essources = null;
    
    /** list of ESModules */
    private ArrayList<ESModuleInstance> esmodules = null;
    
    /** list of Services */
    private ArrayList<ServiceInstance>  services = null;
    
    /** list of Modules */
    private ArrayList<ModuleInstance>   modules = null;
    
    /** list of Paths */
    private ArrayList<Path>             paths = null;
    
    /** list of Sequences */
    private ArrayList<Sequence>         sequences = null;
    
    /** list of streams */
    private ArrayList<Stream>           streams = null;

    /** list of primary datasets */
    private ArrayList<PrimaryDataset>   datasets = null;


    //
    // construction
    //
    
    /** empty constructor */
    public Configuration()
    {
	psets         = new ArrayList<PSetParameter>();
	edsources     = new ArrayList<EDSourceInstance>();
	essources     = new ArrayList<ESSourceInstance>();
	esmodules     = new ArrayList<ESModuleInstance>();
	services      = new ArrayList<ServiceInstance>();
	modules       = new ArrayList<ModuleInstance>();
	paths         = new ArrayList<Path>();
	sequences     = new ArrayList<Sequence>();
	streams       = new ArrayList<Stream>();
	datasets      = new ArrayList<PrimaryDataset>();
    }
    
    /** standard constructor */
    public Configuration(ConfigInfo configInfo,SoftwareRelease release)
    {
	psets         = new ArrayList<PSetParameter>();
	edsources     = new ArrayList<EDSourceInstance>();
	essources     = new ArrayList<ESSourceInstance>();
	esmodules     = new ArrayList<ESModuleInstance>();
	services      = new ArrayList<ServiceInstance>();
	modules       = new ArrayList<ModuleInstance>();
	paths         = new ArrayList<Path>();
	sequences     = new ArrayList<Sequence>();
	streams       = new ArrayList<Stream>();
	datasets      = new ArrayList<PrimaryDataset>();
	
	initialize(configInfo,release);
    }
    
    
    //
    // public member functions
    //

    /** new configuration*/
    public void initialize(ConfigInfo configInfo,SoftwareRelease release)
    {
	this.configInfo  = configInfo;
	this.release     = release;
	
	setHasChanged(false);

	psets.clear();
	edsources.clear();
	essources.clear();
	services.clear();
	modules.clear();
	paths.clear();
	sequences.clear();
	streams.clear();
	datasets.clear();
    }

    /** reset configuration */
    public void reset()
    { 
	configInfo = null;
	release    = null;
	setHasChanged(false);
	
	psets.clear();
	edsources.clear();
	essources.clear();
	services.clear();
	modules.clear();
	paths.clear();
	sequences.clear();
	streams.clear();
	datasets.clear();
    }
    
    /** set the configuration info */
    public void setConfigInfo(ConfigInfo configInfo)
    {
	if (!configInfo.releaseTag().equals(releaseTag())) {
	    System.err.println("Configuration.setConfigInfo ERROR: "+
			       "releaseTag mismatch (" +
			       releaseTag() + " / " +
			       configInfo.releaseTag() + ")");
	}
	this.configInfo = configInfo;
    }
    
    /** overlaod toString() */
    public String toString()
    {
	String result=new String();
	if (configInfo==null) return result;
	if (parentDir()!=null) result += parentDir().name();
	if (result.length()!=1) result += "/";
	result += name() + "/V" + version();
	return result;
    }
    
    /** number of components of a certain type */
    public int componentCount(Class<?> c)
    {
	if      (c == PSetParameter.class)    return psetCount();
	else if (c == EDSourceInstance.class) return edsourceCount(); 
	else if (c == ESSourceInstance.class) return essourceCount(); 
	else if (c == ESModuleInstance.class) return esmoduleCount(); 
	else if (c == ServiceInstance.class)  return serviceCount(); 
	else if (c == Path.class)             return pathCount(); 
	else if (c == Sequence.class)         return sequenceCount();
	else if (c == ModuleInstance.class)   return moduleCount();
	else if (c == Stream.class)           return streamCount();
	else if (c == PrimaryDataset.class)   return datasetCount();
	System.err.println("ERROR: unknwon class " + c.getName());
	return 0;
    }

    /** isEmpty() */
    public boolean isEmpty()
    {
	return (name().length()==0&&psets.isEmpty()&&
		edsources.isEmpty()&&essources.isEmpty()&&
		services.isEmpty()&&modules.isEmpty()&&
		paths.isEmpty()&&sequences.isEmpty()&&
		streams.isEmpty()&&datasets.isEmpty());
    }

    /** check if configuration and all its versions are locked */
    public boolean isLocked()
    {
	return (configInfo!=null) ? configInfo.isLocked() : false;
    }

    /** check by which user the configuration and all its versions are locked */
    public String lockedByUser()
    {
	return (configInfo!=null) ? configInfo.lockedByUser() : new String();
    }

    /** database identifier */
    public int dbId()
    {
	return (configInfo!=null) ? configInfo.dbId() : -1;
    }

    /** get configuration name */
    public String name()
    {
	return (configInfo!=null) ? configInfo.name() : "";
    }
    
    /** get parent directory */
    public Directory parentDir()
    {
	return (configInfo!=null) ? configInfo.parentDir() : null;
    }

    /** get parent directory database id */
    public int parentDirId() 
    {
	return (parentDir()!=null) ? parentDir().dbId() : 0;
    }
    
    /** get configuration version */
    public int version()
    {
	return (configInfo!=null) ? configInfo.version() : 0;
    }
    
    /** next version */
    public int nextVersion()
    {
	return (configInfo!=null) ? configInfo.nextVersion() : 0;
    }
    
    /** add the next version */
    public void addNextVersion(int versionId,
			       String created,String creator,
			       String releaseTag,String processName,
			       String comment)
    {
	configInfo.addVersion(versionId,nextVersion(),created,creator,
			      releaseTag,processName,comment);
	configInfo.setVersionIndex(0);
    }
    
    /** get configuration date of creation as a string */
    public String created()
    {
	return (configInfo!=null) ? configInfo.created() : "";
    }
    
    /** get configuration creator */
    public String creator()
    {
	return (configInfo!=null) ? configInfo.creator() : "";
    }
    
    /** get release tag this configuration is associated with */
    public String releaseTag()
    {
	return (configInfo!=null) ? configInfo.releaseTag() : "";
    }

    /** get the process name */
    public String processName()
    {
	return (configInfo!=null) ? configInfo.processName() : "";
    }
    
    /** get the comment */
    public String comment()
    {
	return (configInfo!=null) ? configInfo.comment() : "";
    }
    
    /** get the software release */
    public SoftwareRelease release() { return release; }

    /** indicate if configuration must be saved */
    public boolean hasChanged()
    {
	if (hasChanged) return true;
	for (EDSourceInstance eds : edsources) if (eds.hasChanged()) return true;
	for (ESSourceInstance ess : essources) if (ess.hasChanged()) return true;
	for (ESModuleInstance esm : esmodules) if (esm.hasChanged()) return true;
	for (ServiceInstance  svc : services)  if (svc.hasChanged()) return true;
	for (Path             pth : paths)     if (pth.hasChanged()) return true;
	for (Sequence         seq : sequences) if (seq.hasChanged()) return true;
	for (Stream           stm : streams)   if (stm.hasChanged()) return true;
	for (PrimaryDataset   dat : datasets)  if (dat.hasChanged()) return true;
	return false;
    }
    
    /** set the 'hasChanged' flag */
    public void setHasChanged(boolean hasChanged) { this.hasChanged = hasChanged; }
    
    /** check if a qualifier is unique */
    public boolean isUniqueQualifier(String qualifier)
    {
	if (qualifier.length()==0) return false;
	for (ESSourceInstance ess : essources)
	    if (ess.name().equals(qualifier)) return false;
	for (ESModuleInstance esm : esmodules)
	    if (esm.name().equals(qualifier)) return false;
	for (ModuleInstance m : modules)
	    if (m.name().equals(qualifier)) return false;
	for (Path p : paths)
	    if (p.name().equals(qualifier)) return false;
	for (Sequence s : sequences)
	    if (s.name().equals(qualifier)) return false;
	return true;
    }
    
    /** check if the reference container has a unique qualifier */
    public boolean hasUniqueQualifier(Referencable referencable)
    {
	if (referencable.name().length()==0) return false;
	for (ESSourceInstance ess : essources)
	    if (ess.name().equals(referencable.name())) return false;
	for (ESModuleInstance esm : esmodules)
	    if (esm.name().equals(referencable.name())) return false;
	for (ModuleInstance m : modules) {
	    if (m==referencable) continue;
	    if (m.name().equals(referencable.name())) return false;
	}
	for (Path p : paths) {
	    if (p==referencable) continue;
	    if (p.name().equals(referencable.name())) return false;
	}
	for (Sequence s : sequences) {
	    if (s==referencable) continue;
	    if (s.name().equals(referencable.name())) return false;
	}
	return true;
    }

    /** check if all entries of a reference container are unique */
    public boolean hasUniqueEntries(ReferenceContainer container)
    {
	for (int i=0;i<container.entryCount();i++) {
	    Reference entry = container.entry(i);
	    if (entry.parent() instanceof ReferenceContainer) {
		ReferenceContainer c = (ReferenceContainer)entry.parent();
		if (!hasUniqueQualifier(c)) return false;
		if (!hasUniqueEntries(c)) return false;
	    }
	    else if (!isUniqueQualifier(entry.name())) return false;
	}
	return true;
    }

    
    /** number of empty containers (paths / sequences) */
    public int emptyContainerCount()
    {
	int result = 0;
	Iterator<Path> itP = paths.iterator();
	while (itP.hasNext()) {
	    Path p = itP.next();
	    if (p.entryCount()==0) result++;
	}
	Iterator<Sequence> itS = sequences.iterator();
	while (itS.hasNext()) {
	    Sequence s = itS.next();
	    if (s.entryCount()==0) result++;
	}
	return result;
    }
    

    //
    // unset tracked parameter counts
    //

    /** total number of unset tracked parameters */
    public int unsetTrackedParameterCount()
    {
	int result = 0;
	result += unsetTrackedPSetParameterCount();
	result += unsetTrackedEDSourceParameterCount();
	result += unsetTrackedESSourceParameterCount();
	result += unsetTrackedESModuleParameterCount();
	result += unsetTrackedServiceParameterCount();
	result += unsetTrackedModuleParameterCount();
	return result;
    }

    /** number of unsert tracked global pset parameters */
    public int unsetTrackedPSetParameterCount()
    {
	int result = 0;
	for (PSetParameter pset : psets)
	    result += pset.unsetTrackedParameterCount();
	return result;
    }
    
    /** number of unsert tracked edsource parameters */
    public int unsetTrackedEDSourceParameterCount()
    {
	int result = 0;
	for (EDSourceInstance eds : edsources)
	    result+=eds.unsetTrackedParameterCount();
	return result;
    }

    /** number of unsert tracked essource parameters */
    public int unsetTrackedESSourceParameterCount()
    {
	int result = 0;
	for (ESSourceInstance ess : essources)
	    result+=ess.unsetTrackedParameterCount();
	return result;
    }

    /** number of unsert tracked esmodule parameters */
    public int unsetTrackedESModuleParameterCount()
    {
	int result = 0;
	for (ESModuleInstance esm : esmodules)
	    result+=esm.unsetTrackedParameterCount();
	return result;
    }

    /** number of unsert tracked service parameters */
    public int unsetTrackedServiceParameterCount()
    {
	int result = 0;
	for (ServiceInstance svc : services)
	    result+=svc.unsetTrackedParameterCount();
	return result;
    }

    /** number of unsert tracked module parameters */
    public int unsetTrackedModuleParameterCount()
    {
	int result = 0;
	for (ModuleInstance mod : modules)
	    result+=mod.unsetTrackedParameterCount();
	return result;
    }

    /** number of paths unassigned to any stream */
    public int pathNotAssignedToStreamCount()
    {
	int result = 0;
	if (streams.size()==0) return result;
	for (Path p : paths) {
	    if (p.isEndPath()) continue;
	    if (p.streamCount()==0) result++;
	}
	return result;
    }

    /** number of paths unassigned to any primary dataset */
    public int pathNotAssignedToDatasetCount()
    {
	int result = 0;
	if (datasets.size()==0) return result;
	for (Path p : paths) {
	    if (p.isEndPath()) continue;
	    if (p.datasetCount()==0) result++;
	}
	return result;
    }


    //
    // PSets
    //
    
    /**  number of global PSets */
    public int psetCount() { return psets.size(); }

    /** get i-th global PSet */
    public PSetParameter pset(int i) { return psets.get(i); }

    /** get global pset by name */
    public PSetParameter pset(String name)
    {
	for (PSetParameter pset : psets)
	    if (pset.name().equals(name)) return pset;
	return null;
    }

    /** index of a certain global PSet */
    public int indexOfPSet(PSetParameter pset)
    {
	return psets.indexOf(pset);
    }

    /** retrieve pset iterator */
    public Iterator<PSetParameter> psetIterator() { return psets.iterator(); }
    
    /** insert global pset at i-th position */
    public void insertPSet(PSetParameter pset)
    {
	psets.add(pset);
	hasChanged = true;
    }
    
    /** remove a global PSet */
    public void removePSet(PSetParameter pset)
    {
	psets.remove(pset);
	hasChanged = true;
    }
    
    /** sort global PSets*/
    public void sortPSets() { Collections.sort(psets); hasChanged=true; }
    

    //
    // EDSources 
    //
    
    /**  number of EDSources */
    public int edsourceCount() { return edsources.size(); }

    /** get i-th EDSource */
    public EDSourceInstance edsource(int i) { return edsources.get(i); }

    /** get EDSource by name */
    public EDSourceInstance edsource(String name)
    {
	for (EDSourceInstance eds : edsources)
	    if (eds.name().equals(name)) return eds;
	return null;
    }

    /** index of a certain EDSource */
    public int indexOfEDSource(EDSourceInstance edsource)
    {
	return edsources.indexOf(edsource);
    }
    
    /** retrieve edsource iterator */
    public Iterator<EDSourceInstance> edsourceIterator() { return edsources.iterator(); }

    /** insert EDSource at i-th position */
    public EDSourceInstance insertEDSource(String templateName)
    {
	if (edsourceCount()>0) return null;

	EDSourceTemplate template =
	    (EDSourceTemplate)release.edsourceTemplate(templateName);
	if (template==null) {
	    System.err.println("insertEDSource ERROR: unknown template '" +
			       templateName+"'!");
	    return null;
	}

	EDSourceInstance instance = null;	
	try {
	    instance = (EDSourceInstance)template.instance();
	    edsources.add(instance);
	    instance.setConfiguration(this);
	    hasChanged = true;
	}
	catch (Exception e) {
	    System.out.println(e.getMessage());
	}
	return instance;
    }
    
    /** remove a EDSource */
    public void removeEDSource(EDSourceInstance edsource)
    {
	edsource.remove();
	int index = edsources.indexOf(edsource);
	edsources.remove(index);
	hasChanged = true;
    }
    
    /** sort  EDSources */
    public void sortEDSources() { Collections.sort(edsources); hasChanged=true; }
    

    //
    // ESSources
    //
    
    /**  number of ESSources */
    public int essourceCount() { return essources.size(); }
    
    /** get i-th ESSource */
    public ESSourceInstance essource(int i) { return essources.get(i); }

    /** get ESSource by name */
    public ESSourceInstance essource(String name)
    {
	for (ESSourceInstance ess : essources)
	    if (ess.name().equals(name)) return ess;
	return null;
    }

    /** index of a certain ESSource */
    public int indexOfESSource(ESSourceInstance essource)
    {
	return essources.indexOf(essource);
    }
    
    /** retrieve essource iterator */
    public Iterator<ESSourceInstance> essourceIterator() { return essources.iterator(); }

    /** insert ESSource at i=th position */
    public ESSourceInstance insertESSource(int i,
					   String templateName,
					   String instanceName)
    {
	ESSourceTemplate template =
	    (ESSourceTemplate)release.essourceTemplate(templateName);
	if (template==null) {
	    System.err.println("insertESSource ERROR: unknown template '"+
			       templateName+"'!");
	    return null;
	}
	
	ESSourceInstance instance = null;	
	try {
	    instance = (ESSourceInstance)template.instance(instanceName);
	    essources.add(i,instance);
	    instance.setConfiguration(this);
	    hasChanged = true;
	}
	catch (Exception e) {
	    System.out.println(e.getMessage());
	}
	return instance;
    }
    
    /** remove a ESSource */
    public void removeESSource(ESSourceInstance essource)
    {
	essource.remove();
	int index = essources.indexOf(essource);
	essources.remove(index);
	hasChanged = true;
    }

    /** sort  ESSources */
    public void sortESSources() { Collections.sort(essources); hasChanged=true; }
    
    
    //
    // ESModules
    //
    
    /**  number of ESModules */
    public int esmoduleCount() { return esmodules.size(); }
    
    /** get i-th ESModule */
    public ESModuleInstance esmodule(int i) { return esmodules.get(i); }
    
    /** get ESModule by name */
    public ESModuleInstance esmodule(String name)
    {
	for (ESModuleInstance esm : esmodules)
	    if (esm.name().equals(name)) return esm;
	return null;
    }

    /** index of a certain ESSource */
    public int indexOfESModule(ESModuleInstance esmodule)
    {
	return esmodules.indexOf(esmodule);
    }
   
    /** retrieve esmodule iterator */
    public Iterator<ESModuleInstance> esmoduleIterator() { return esmodules.iterator(); }


    /** insert ESModule at i-th position */
    public ESModuleInstance insertESModule(int i,
					   String templateName,
					   String instanceName)
    {
	ESModuleTemplate template =
	    (ESModuleTemplate)release.esmoduleTemplate(templateName);
	if (template==null) {
	    System.err.println("insertESModule ERROR: unknown template '" +
			       templateName+"'!");
	    return null;
	}

	ESModuleInstance instance = null;
	try {
	    instance = (ESModuleInstance)template.instance(instanceName);
	    esmodules.add(i,instance);
	    instance.setConfiguration(this);
	    hasChanged = true;
	}
	catch (Exception e) {
	    System.out.println(e.getMessage());
	}
	return instance;
    }
    
    /** remove a ESModule */
    public void removeESModule(ESModuleInstance esmodule)
    {
	esmodule.remove();
	int index = esmodules.indexOf(esmodule);
	esmodules.remove(index);
	hasChanged = true;
    }
    
    /** sort  ESModules */
    public void sortESModules() { Collections.sort(esmodules); hasChanged=true; }
    
    
    //
    // Services
    //
    
    /**  number of Services */
    public int serviceCount() { return services.size(); }

    /** get i-th Service */
    public ServiceInstance service(int i) { return services.get(i); }

    /** get Service by name */
    public ServiceInstance service(String name)
    {
	for (ServiceInstance svc : services)
	    if (svc.name().equals(name)) return svc;
	return null;
    }
    
    /** index of a certain Service */
    public int indexOfService(ServiceInstance service)
    {
	return services.indexOf(service);
    }

    /** retrieve service iterator */
    public Iterator<ServiceInstance> serviceIterator() { return services.iterator(); }
    
    /** insert Service at i=th position */
    public ServiceInstance insertService(int i,String templateName)
    {
	ServiceTemplate template =
	    (ServiceTemplate)release.serviceTemplate(templateName);
	if (template==null) {
	    System.err.println("insertService ERROR: unknown template '" +
			       templateName+"'!");
	    return null;
	}

	ServiceInstance instance = null;
	try {
	    instance = (ServiceInstance)template.instance();
	    services.add(i,instance);
	    instance.setConfiguration(this);
	    hasChanged = true;
	}
	catch (Exception e) {
	    System.out.println(e.getMessage());
	}
	return instance;
    }
    
    /** remove a Service */
    public void removeService(ServiceInstance service)
    {
	service.remove();
	int index = services.indexOf(service);
	services.remove(index);
	hasChanged = true;
    }
    
    /** sort services */
    public void sortServices() { Collections.sort(services); hasChanged=true; }
    
    
    //
    // Modules 
    //
    
    /**  number of Modules */
    public int moduleCount() { return modules.size(); }

    /** get i-th Module */
    public ModuleInstance module(int i) { return modules.get(i); }
    
    /** get Module by name */
    public ModuleInstance module(String moduleName)
    {
	for (ModuleInstance m : modules) if (m.name().equals(moduleName)) return m;
	return null;
    }
    
    /** index of a certain Module */
    public int indexOfModule(ModuleInstance module)
    {
	return modules.indexOf(module);
    }
    
    /** retrieve module iterator */
    public Iterator<ModuleInstance> moduleIterator() { return modules.iterator(); }
    
    /** insert a module */
    public ModuleInstance insertModule(String templateName,String instanceName)
    {
	ModuleTemplate template =
	    (ModuleTemplate)release.moduleTemplate(templateName);
	if (template == null) {
	    System.err.println("insertModule ERROR: unknown template '" +
			       templateName+"' (instanceName="+instanceName+")!");
	    return null;
	}

	ModuleInstance instance = null;
	try {
	    instance = (ModuleInstance)template.instance(instanceName);
	    if (instance.referenceCount()==0) {
		modules.add(instance);
		instance.setConfiguration(this);
		hasChanged = true;
	    }
	}
	catch (Exception e) {
	    System.out.println(e.getMessage());
	}
	return instance;
    }

    /** insert a pre-existing module */
    public boolean insertModule(int i,ModuleInstance module)
    {
	if (modules.indexOf(module)<0&&module.referenceCount()==0) {
	    modules.add(i,module);
	    module.setConfiguration(this);
	    hasChanged = true;
	    return true;
	}
	return false;
    }
    
    /** remove a module reference */
    public void removeModuleReference(ModuleReference module)
    {
	ModuleInstance instance = (ModuleInstance)module.parent();
	module.remove();
	if (instance.referenceCount()==0) {
	    int index = modules.indexOf(instance);
	    modules.remove(index);
	}
	hasChanged = true;
    }
    
    
    /** insert ModuleReference at i-th position into a path/sequence */
    public ModuleReference insertModuleReference(ReferenceContainer container,
						 int                i,
						 ModuleInstance     instance)
    {
	ModuleReference reference =
	    (ModuleReference)instance.createReference(container,i);
	hasChanged = true;
	return reference;
    }
    
    /** insert ModuleReference at i-th position into a path/sequence */
    public ModuleReference insertModuleReference(ReferenceContainer container,
						 int                i,
						 String             templateName,
						 String             instanceName)
    {
	ModuleInstance instance = insertModule(templateName,instanceName);
	return (instance!=null) ?
	    insertModuleReference(container,i,instance) : null;
    }    
    
    /** sort  Modules */
    public void sortModules() { Collections.sort(modules); }

    
    //
    // Paths
    //

    /** number of Paths */
    public int pathCount() { return paths.size(); }
    
    /** get i-th Path */
    public Path path(int i) { return paths.get(i); }

    /** get Path by name*/
    public Path path(String pathName)
    {
	for (Path p : paths)
	    if (p.name().equals(pathName)) return p;
	return null;
    }

    /** index of a certain Path */
    public int indexOfPath(Path path) { return paths.indexOf(path); }
    
    /** retrieve path iterator */
    public Iterator<Path> pathIterator() { return paths.iterator(); }

    /** insert path at i-th position */
    public Path insertPath(int i, String pathName)
    {
	Path path = new Path(pathName);
	paths.add(i,path);
	path.setConfiguration(this);
	hasChanged = true;
	return path;
    }
    
    /** move a path to another position within paths */
    public boolean movePath(Path path,int targetIndex)
    {
	int currentIndex = paths.indexOf(path);
	if (currentIndex<0) return false;
	if (currentIndex==targetIndex) return true;
	if (targetIndex>paths.size()) return false;
	if (currentIndex<targetIndex) targetIndex--;
	paths.remove(currentIndex);
	paths.add(targetIndex,path);
	hasChanged = true;
	return true;
    }
    
    /** get the sequence number of a certain path */
    public int pathSequenceNb(Path path) { return paths.indexOf(path); }
    
    /** remove a path */
    public void removePath(Path path)
    {
	while (path.referenceCount()>0) {
	    PathReference reference = (PathReference)path.reference(0);
	    reference.remove();
	}
	
	// remove all entries of this path
	while (path.entryCount()>0) {
	    Reference reference = path.entry(0);
	    reference.remove();
	    if (reference instanceof ModuleReference) {
		ModuleReference module   = (ModuleReference)reference;
		ModuleInstance  instance = (ModuleInstance)module.parent();
		if (instance.referenceCount()==0) {
		    int index = modules.indexOf(instance);
		    modules.remove(index);
		}
	    }
	}
	
	// remove this paths from all parent primary datasets
	Iterator<PrimaryDataset> itPD = path.datasetIterator();
	while (itPD.hasNext()) {
	    PrimaryDataset pd = itPD.next();
	    itPD.remove();
	    pd.removePath(path);
	}
	
	int index = paths.indexOf(path);
	paths.remove(index);
	hasChanged = true;
    }
    
    /** insert a path reference into another path/sequence */
    public PathReference insertPathReference(ReferenceContainer parentPath,
					     int i,Path path)
    {
	PathReference reference = (PathReference)path.createReference(parentPath,i);
	hasChanged = true;
	return reference;
    }
    
    /** sort Paths */
    public void sortPaths() { Collections.sort(paths); hasChanged=true; }

    
    //
    // Sequences
    //

    /** number of Sequences */
    public int sequenceCount() { return sequences.size(); }
    
    /** get i-th Sequence */
    public Sequence sequence(int i) { return sequences.get(i); }

    /** get Sequence by name*/
    public Sequence sequence(String sequenceName)
    {
	for (Sequence s : sequences)
	    if (s.name().equals(sequenceName)) return s;
	return null;
    }
    
    /** index of a certain Sequence */
    public int indexOfSequence(Sequence sequence)
    {
	return sequences.indexOf(sequence);
    }
    
    /** retrieve sequence iterator */
    public Iterator<Sequence> sequenceIterator() { return sequences.iterator(); }
    
    /** retrieve sequence iterator */
    public Iterator<Sequence> orderedSequenceIterator() { return sequenceIterator(); }
    
    /** insert sequence */
    public Sequence insertSequence(int i,String sequenceName)
    {
	Sequence sequence = new Sequence(sequenceName);
	sequences.add(i,sequence);
	sequence.setConfiguration(this);
	hasChanged = true;
	return sequence;
    }
    
    /** move a sequence to another position within sequences */
    public boolean moveSequence(Sequence sequence,int targetIndex)
    {
	int currentIndex = sequences.indexOf(sequence);
	if (currentIndex<0) return false;
	if (currentIndex==targetIndex) return true;
	if (targetIndex>sequences.size()) return false;
	if (currentIndex<targetIndex) targetIndex--;
	sequences.remove(currentIndex);
	sequences.add(targetIndex,sequence);
	hasChanged = true;
	return true;
    }
    
    /** remove a sequence */
    public void removeSequence(Sequence sequence)
    {
	while (sequence.referenceCount()>0) {
	    SequenceReference reference = (SequenceReference)sequence.reference(0);
	    reference.remove();
	}
	
	// remove all modules from this sequence
	while (sequence.entryCount()>0) {
	    Reference reference = sequence.entry(0);
	    reference.remove();
	    if (reference instanceof ModuleReference) {
		ModuleReference module   = (ModuleReference)reference;
		ModuleInstance  instance = (ModuleInstance)module.parent();
		if (instance.referenceCount()==0) {
		    int index = modules.indexOf(instance);
		    modules.remove(index);
		}
	    }
	}
	
	int index = sequences.indexOf(sequence);
	sequences.remove(index);
	hasChanged = true;
    }
    
    /** insert a sequence reference into another path */
    public SequenceReference insertSequenceReference(ReferenceContainer parent,int i,
						     Sequence sequence)
    {
	SequenceReference reference =
	    (SequenceReference)sequence.createReference(parent,i);
	hasChanged = true;
	return reference;
    }

    /** sort Sequences */
    public void sortSequences() { Collections.sort(sequences); hasChanged=true; }
    

    //
    // Streams
    //
    
    /** number of streams */
    public int streamCount() { return streams.size(); }
    
    /** retrieve i-th stream */
    public Stream stream(int i) { return streams.get(i); }

    /** retrieve stream by label */
    public Stream stream(String streamLabel)
    {
	for (Stream s : streams)
	    if (s.label().equals(streamLabel)) return s;
	return null;
    }

    /** index of a certain stream */
    public int indexOfStream(Stream stream) { return streams.indexOf(stream); }

    /** retrieve stream iterator */
    public Iterator<Stream> streamIterator() { return streams.iterator(); }
    
    /** insert a new stream */
    public Stream insertStream(String streamLabel)
    {
	Stream stream = new Stream(streamLabel);
	streams.add(stream);
	hasChanged=true;
	return stream;
    }

    /** remove a stream */
    public void removeStream(Stream stream)
    {
	if (streams.indexOf(stream)>=0) {
	    Iterator<PrimaryDataset> itD = stream.datasetIterator();
	    while (itD.hasNext()) itD.next().removeFromStream(stream);
	    streams.remove(stream);
	    hasChanged = true;
	}
    }

    /** sort Streams */
    public void sortStreams() { Collections.sort(streams); hasChanged=true;}
    
    
    //
    // Primary Datasets
    //
    
    /** number of primary datasets */
    public int datasetCount() { return datasets.size(); }
    
    /** retrieve i-th primary dataset */
    public PrimaryDataset dataset(int i) { return datasets.get(i); }

    /** retrieve primary dataset by label */
    public PrimaryDataset dataset(String datasetLabel)
    {
	for (PrimaryDataset pd : datasets)
	    if (pd.label().equals(datasetLabel)) return pd;
	return null;
    }
    
    /** index of a certain primary dataset */
    public int indexOfDataset(PrimaryDataset dataset)
    {
	return datasets.indexOf(dataset);
    }

    /** retrieve primary dataset iterator */
    public Iterator<PrimaryDataset> datasetIterator()
    {
	return datasets.iterator();
    }
    
    /** insert a new primary dataset */
    public PrimaryDataset insertDataset(String datasetLabel)
    {
	PrimaryDataset dataset = new PrimaryDataset(datasetLabel);
	datasets.add(dataset);
	hasChanged=true;
	return dataset;
    }

    /** remove a dataset */
    public void removeDataset(PrimaryDataset dataset)
    {
	if (datasets.indexOf(dataset)>=0) {
	    while (dataset.pathCount()>0) dataset.removePath(dataset.path(0));
	    datasets.remove(dataset);
	    hasChanged = true;
	}
    }
    
    /** sort primary datasets */
    public void sortDatasets() { Collections.sort(datasets); hasChanged=true; }
    
}
