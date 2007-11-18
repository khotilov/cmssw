package confdb.data;

import java.util.ArrayList;
import java.util.Iterator;


/**
 * ConfigurationModifier
 * ---------------------
 * @author Philipp Schieferdecker
 *
 * retrieve only selected parts of a configuraiton,
 * after applying configurable filter decisions.
 */
public class ConfigurationModifier implements IConfiguration
{
    //
    // member data
    //

    /** the master configuration */
    private IConfiguration master = null;

    /** flag indicating if the filter was applied */
    private boolean isModified = false;
    
    /** filtered components */
    private ArrayList<PSetParameter>    psets    =new ArrayList<PSetParameter>();
    private ArrayList<EDSourceInstance> edsources=new ArrayList<EDSourceInstance>();
    private ArrayList<ESSourceInstance> essources=new ArrayList<ESSourceInstance>();
    private ArrayList<ESModuleInstance> esmodules=new ArrayList<ESModuleInstance>();
    private ArrayList<ServiceInstance>  services =new ArrayList<ServiceInstance>();
    private ArrayList<ModuleInstance>   modules  =new ArrayList<ModuleInstance>();
    private ArrayList<Path>             paths    =new ArrayList<Path>();
    private ArrayList<Sequence>         sequences=new ArrayList<Sequence>();
    private ArrayList<Stream>           streams  =new ArrayList<Stream>();

    /** internal instructions */
    private ModifierInstructions modifications = new ModifierInstructions();
    
    //
    // construction
    //

    /** standard constructor */
    public ConfigurationModifier(IConfiguration config)
    {
	this.master = config;
    }

    //
    // member functions
    //

    /** replace the current EDSource with a PoolSource */
    public void insertPoolSource(String fileNames)
    {
	modifications.insertPoolSource(fileNames);
    }

    /** replace the current EDSource with a DaqSource */
    public void insertDaqSource()
    {
	if (edsourceCount()==0||!edsource(0).template().name().equals("DaqSource"))
	    modifications.insertDaqSource();
    }
    
    /** replace current OutputModules with ShmStreamConsumer */
    public void insertShmStreamConsumer() { modifications.insertShmStreamConsumer();}

    /** replace current OutputModules with PoolOutputModules */
    public void insertPoolOutputModule(String fileName)
    {
	modifications.insertPoolOutputModule(fileName);
    }

    /** apply modifications based on internal modifications */
    public void modify()
    {
	modify(modifications);
    }

    /** apply modifications */
    public void modify(ModifierInstructions modifications)
    {
	if (modifications==null) return;
	
	psets.clear();
	edsources.clear();
	essources.clear();
	esmodules.clear();
	services.clear();
	modules.clear();
	paths.clear();
	sequences.clear();
	streams.clear();
	
	if (!modifications.resolve(master)) return;
	
	if (!modifications.doFilterAllPSets()) {
	    Iterator<PSetParameter> it=master.psetIterator();
	    while (it.hasNext()) {
		PSetParameter pset = it.next();
		if (!modifications.isInBlackList(pset))
		    psets.add(pset);
	    }
	}
	
	if (!modifications.doFilterAllEDSources()) {
	    Iterator<EDSourceInstance> it=master.edsourceIterator();
	    while (it.hasNext()) {
		EDSourceInstance edsource = it.next();
		if (!modifications.isInBlackList(edsource))
		    edsources.add(edsource);
	    }
	}
	if (modifications.doInsertEDSource())
	    edsources.add(modifications.edsourceToBeAdded());
	
	if (!modifications.doFilterAllESSources()) {
	    Iterator<ESSourceInstance> it=master.essourceIterator();
	    while (it.hasNext()) {
		ESSourceInstance essource = it.next();
		if (!modifications.isInBlackList(essource))
		    essources.add(essource);
	    }
	}
	
	if (!modifications.doFilterAllESModules()) {
	    Iterator<ESModuleInstance> it=master.esmoduleIterator();
	    while (it.hasNext()) {
		ESModuleInstance esmodule = it.next();
		if (!modifications.isInBlackList(esmodule))
		    esmodules.add(esmodule);
	    }
	}
	
	if (!modifications.doFilterAllServices()) {
	    Iterator<ServiceInstance> it=master.serviceIterator();
	    while (it.hasNext()) {
		ServiceInstance service = it.next();
		if (!modifications.isInBlackList(service))
		    services.add(service);
	    }
	}
	
	boolean hasOutputModule = false;
	
	if (!modifications.doFilterAllPaths()) {
	    Iterator<Path> itP = master.pathIterator();
	    while (itP.hasNext()) {
		Path path = itP.next();
		if (!modifications.isInBlackList(path)) {
		    
		    if (path.hasOutputModule()) hasOutputModule = true;
		    
		    if (path.hasOutputModule()&&
			modifications.doInsertOutputModule()) {

			Path copy = new Path(path.name());
			Iterator<Reference> it = path.entryIterator();
			while (it.hasNext()) {
			    Reference entry = it.next();
			    ModuleInstance outputModule = null;
			    if (entry instanceof ModuleReference) {
				ModuleReference ref  = (ModuleReference)entry;
				ModuleInstance  inst = (ModuleInstance)ref.parent();
				if (inst.template().type().equals("OutputModule"))
				    outputModule = inst;
			    }
			    if (outputModule!=null) {
				ModuleInstance outputI = 
				    modifications.outputModuleToBeAdded(entry
								       .name());
				outputI
				    .updateParameter("SelectEvents","PSet",
						     outputModule
						     .parameter("SelectEvents",
								"PSet")
						     .valueAsString());
				outputI
				    .updateParameter("outputCommands","vstring",
						     outputModule
						     .parameter("outputCommands",
								"vstring")
						     .valueAsString());
				
				outputI.createReference(copy,copy.entryCount());
			    }
			    else
				copy.insertEntry(copy.entryCount(),entry);
			}
			path = copy;
		    }

		    paths.add(path);
		    Iterator<Sequence> itS = path.sequenceIterator();
		    while (itS.hasNext()) {
			Sequence sequence = itS.next();
			if (sequences.indexOf(sequence)<0)
			    sequences.add(sequence);
		    }
		    Iterator<ModuleInstance> itM = path.moduleIterator();
		    while (itM.hasNext()) {
			ModuleInstance module = itM.next();
			if (modules.indexOf(module)<0)
			    modules.add(module);
		    }
		}
	    }
	}
	
	if (!hasOutputModule&&modifications.doInsertOutputModule()) {
	    Path out = new Path("output");
	    ModuleInstance outputI = modifications.outputModuleToBeAdded("out");
	    outputI.createReference(out,0);
	    paths.add(out);
	}

	Iterator<String> itS = modifications.requestedSequenceIterator();
	while (itS.hasNext()) {
	    String   sequenceName = itS.next();
	    Sequence sequence     = master.sequence(sequenceName);
	    if (sequence!=null&&sequences.indexOf(sequence)<0)
		sequences.add(sequence);
	}
	


	Iterator<String> itM = modifications.requestedModuleIterator();
	while (itM.hasNext()) {
	    String         moduleLabel = itM.next();
	    ModuleInstance module      = master.module(moduleLabel);
	    if (module!=null&&modules.indexOf(module)<0)
		modules.add(module);
	}
	
	isModified = true;

	System.out.println("sequenceCount = " + sequenceCount());
	System.out.println("moduleCount   = " + moduleCount());
    }

    /** reset all modifications */
    public void reset()
    {
	modifications = new ModifierInstructions();
	
	psets.clear();
	edsources.clear();
	essources.clear();
	esmodules.clear();
	services.clear();
	modules.clear();
	paths.clear();
	sequences.clear();
	streams.clear();

	isModified = false;
    }
    
    
    //
    // IConfiguration interface
    //
    
    /** name of the configuration */
    public String name() { return master.name(); }

    /** parent directory of the configuration */
    public Directory parentDir() { return master.parentDir(); }

    /** version of the configuration */
    public int version() { return master.version(); }
    
    /** get configuration date of creation as a string */
    public String created() { return master.created(); }
	
    /** get configuration creator */
    public String creator() { return master.creator(); }
    
    /** release tag */
    public String releaseTag() { return master.releaseTag(); }

    /** process name */
    public String processName() { return master.processName(); }

    /** comment */
    public String comment() { return master.comment(); }
    
    
    /** total number of components of a certain type */
    public int componentCount(Class<?> c)
    {
	return master.componentCount(c);
    }

    /** check if the configuration is empty */
    public boolean isEmpty() { return false; }

    /** check if the provided string is already in use as a label */
    public boolean isUniqueQualifier(String qualifier)
    {
	return master.isUniqueQualifier(qualifier);
    }

    
    /** total number of unset tracked parameters */
    public int unsetTrackedParameterCount() {
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
	if (!isModified) return master.unsetTrackedPSetParameterCount();
	int result = 0;
	for (PSetParameter pset : psets)
	    result += pset.unsetTrackedParameterCount();
	return result;
    }
    
    /** number of unsert tracked edsource parameters */
    public int unsetTrackedEDSourceParameterCount()
    {
	if (!isModified) return master.unsetTrackedEDSourceParameterCount();
	int result = 0;
	for (EDSourceInstance eds : edsources)
	    result+=eds.unsetTrackedParameterCount();
	return result;
    }

    /** number of unsert tracked essource parameters */
    public int unsetTrackedESSourceParameterCount()
    {
	if (!isModified) return master.unsetTrackedESSourceParameterCount();
	int result = 0;
	for (ESSourceInstance ess : essources)
	    result+=ess.unsetTrackedParameterCount();
	return result;
    }

    /** number of unsert tracked esmodule parameters */
    public int unsetTrackedESModuleParameterCount()
    {
	if (!isModified) return master.unsetTrackedESModuleParameterCount();
	int result = 0;
	for (ESModuleInstance esm : esmodules)
	    result+=esm.unsetTrackedParameterCount();
	return result;
    }

    /** number of unsert tracked service parameters */
    public int unsetTrackedServiceParameterCount()
    {
	if (!isModified) return master.unsetTrackedServiceParameterCount();
	int result = 0;
	for (ServiceInstance svc : services)
	    result+=svc.unsetTrackedParameterCount();
	return result;
    }

    /** number of unsert tracked module parameters */
    public int unsetTrackedModuleParameterCount()
    {
	if (!isModified) return master.unsetTrackedModuleParameterCount();
	int result = 0;
	for (ModuleInstance mod : modules)
	    result+=mod.unsetTrackedParameterCount();
	return result;
    }

    /** number of paths unassigned to any stream */
    public int pathNotAssignedToStreamCount()
    {
	if (!isModified) return master.pathNotAssignedToStreamCount();
	int result = 0;
	if (streams.size()==0) return result;
	for (Path p : paths) if (p.streamCount()==0) result++;
	return result;
    }
    
    /**  number of global PSets */
    public int psetCount()
    {
	return (isModified) ? psets.size() : master.psetCount();
    }

    /** get i-th global PSet */
    public PSetParameter pset(int i)
    {
	return (isModified) ? psets.get(i) : master.pset(i);
    }
    
    /** index of a certain global PSet */
    public int indexOfPSet(PSetParameter pset)
    {
	return (isModified) ? psets.indexOf(pset) : master.indexOfPSet(pset);
    }

    /** retrieve pset iterator */
    public Iterator<PSetParameter> psetIterator()
    {
	return (isModified) ? psets.iterator() : master.psetIterator();
    }
    
    
    /**  number of EDSources */
    public int edsourceCount()
    {
	return (isModified) ? edsources.size() : master.edsourceCount();
    }

    /** get i-th EDSource */
    public EDSourceInstance edsource(int i)
    {
	return (isModified) ? edsources.get(i) : master.edsource(i);
    }

    /** index of a certain EDSource */
    public int indexOfEDSource(EDSourceInstance edsource)
    {
	return (isModified) ?
	    edsources.indexOf(edsource) : master.indexOfEDSource(edsource);
    }
	
    /** retrieve edsource iterator */
    public Iterator<EDSourceInstance> edsourceIterator()
    {
	return (isModified) ? edsources.iterator() : master.edsourceIterator();
    }

    
    /**  number of ESSources */
    public int essourceCount()
    {
	return (isModified) ? essources.size() : master.essourceCount();
    }
    
    /** get i-th ESSource */
    public ESSourceInstance essource(int i)
    {
	return (isModified) ? essources.get(i) : master.essource(i);
    }

    /** index of a certain ESSource */
    public int indexOfESSource(ESSourceInstance essource)
    {
	return (isModified) ?
	    essources.indexOf(essource) : master.indexOfESSource(essource);
    }

    /** retrieve essource iterator */
    public Iterator<ESSourceInstance> essourceIterator()
    {
	return (isModified) ? essources.iterator() : master.essourceIterator();
    }
    
    
    /**  number of ESModules */
    public int esmoduleCount()
    {
	return (isModified) ? esmodules.size() : master.esmoduleCount();
    }
    
    /** get i-th ESModule */
    public ESModuleInstance esmodule(int i)
    {
	return (isModified) ? esmodules.get(i) : master.esmodule(i);
    }

    /** index of a certain ESSource */
    public int indexOfESModule(ESModuleInstance esmodule)
    {
	return (isModified) ?
	    esmodules.indexOf(esmodule) : master.indexOfESModule(esmodule);
    }
    
    /** retrieve esmodule iterator */
    public Iterator<ESModuleInstance> esmoduleIterator()
    {
	return (isModified) ? esmodules.iterator() : master.esmoduleIterator();
    }

    
    /**  number of Services */
    public int serviceCount()
    {
	return (isModified) ? services.size() : master.serviceCount();
    }

    /** get i-th Service */
    public ServiceInstance service(int i)
    {
	return (isModified) ? services.get(i) : master.service(i);
    }

    /** index of a certain Service */
    public int indexOfService(ServiceInstance service)
    {
	return (isModified) ?
	    services.indexOf(service) : master.indexOfService(service);
    }
    
    /** retrieve service iterator */
    public Iterator<ServiceInstance> serviceIterator()
    {
	return (isModified) ? services.iterator() : master.serviceIterator();
    }
    
    
    /**  number of Modules */
    public int moduleCount()
    {
	return (isModified) ? modules.size() : master.moduleCount();
    }

    /** get i-th Module */
    public ModuleInstance module(int i)
    {
	return (isModified) ? modules.get(i) : master.module(i);
    }
    
    /** get Module by name */
    public ModuleInstance module(String moduleName)
    {
	return master.module(moduleName);
    }

    /** index of a certain Module */
    public int indexOfModule(ModuleInstance module)
    {
	return (isModified) ?
	    modules.indexOf(module) : master.indexOfModule(module);
    }
    
    /** retrieve module iterator */
    public Iterator<ModuleInstance> moduleIterator()
    {
	return (isModified) ? modules.iterator() : master.moduleIterator();
    }
    

    /** number of Paths */
    public int pathCount()
    {
	return (isModified) ? paths.size() : master.pathCount();
    }
    
    /** get i-th Path */
    public Path path(int i)
    {
	return (isModified) ? paths.get(i) : master.path(i);
    }

    /** geth Path by name */
    public Path path(String pathName)
    {
	return master.path(pathName);
    }

    /** index of a certain Path */
    public int indexOfPath(Path path)
    {
	return (isModified) ? paths.indexOf(path) : master.indexOfPath(path);
    }
    
    /** retrieve path iterator */
    public Iterator<Path> pathIterator()
    {
	return (isModified) ? paths.iterator() : master.pathIterator();
    }

    
    /** number of Sequences */
    public int sequenceCount()
    {
	return (isModified) ? sequences.size() : master.sequenceCount();
    }
    
    /** get i-th Sequence */
    public Sequence sequence(int i)
    {
	return (isModified) ? sequences.get(i) : master.sequence(i);
    }

    /** get Sequence by name */
    public Sequence sequence(String sequenceName)
    {
	return master.sequence(sequenceName);
    }

    /** index of a certain Sequence */
    public int indexOfSequence(Sequence sequence)
    {
	return (isModified) ?
	    sequences.indexOf(sequence) : master.indexOfSequence(sequence);
    }

    /** retrieve sequence iterator */
    public Iterator<Sequence> sequenceIterator()
    {
	return (isModified) ? sequences.iterator() : master.sequenceIterator();
    }
    
    
    /** number of streams */
    public int streamCount()
    {
	return (isModified) ? streams.size() : master.streamCount();
    }
    
    /** retrieve i-th stream */
    public Stream stream(int i)
    {
	return (isModified) ? streams.get(i) : master.stream(i);
    }
    
    /** index of a certain stream */
    public int indexOfStream(Stream stream)
    {
	return (isModified) ?
	    streams.indexOf(stream) : master.indexOfStream(stream);
    }
    
    /** retrieve stream iterator */
    public Iterator<Stream> streamIterator()
    {
	return (isModified) ? streams.iterator() : master.streamIterator();
    }
    
}
