#!/bin/sh
# Pass in name and status
function die { echo $1: status $2 ;  exit $2; }

rm -f ${LOCAL_TMP_DIR}/PoolInputTest.root ${LOCAL_TMP_DIR}/PoolInputOther.root
rm -f ${LOCAL_TMP_DIR}/PoolInputTestCatalog.xml ${LOCAL_TMP_DIR}/PoolInputTestCatalog.xml.BAK
rm -f ${LOCAL_TMP_DIR}/PrePoolInputTest.cfg ${LOCAL_TMP_DIR}/PoolInputTest.cfg

cat > ${LOCAL_TMP_DIR}/PrePoolInputTest.cfg << !
# Configuration file for PrePoolInputTest 
process TEST = {
	path p = {Thing, OtherThing}
	module Thing = ThingProducer {untracked int32 debugLevel = 1}
	module OtherThing = OtherThingProducer {untracked int32 debugLevel = 1}
	module output = PoolOutputModule {
		untracked string fileName = '${LOCAL_TMP_DIR}/PoolInputTest.root'
		untracked string catalog = '${LOCAL_TMP_DIR}/PoolInputTestCatalog.xml'
		untracked string logicalFileName = 'PoolTest.root'
		untracked int32 maxSize = 100000
	}
	source = EmptySource {untracked int32 maxEvents = 2}
        endpath ep = {output}
}
!
cmsRun --parameter-set ${LOCAL_TMP_DIR}/PrePoolInputTest.cfg || die 'Failure using PrePoolInputTest.cfg' $?

cp ${LOCAL_TMP_DIR}/PoolInputTest.root ${LOCAL_TMP_DIR}/PoolInputOther.root

cat > ${LOCAL_TMP_DIR}/PoolInputTest.cfg << !
# Configuration file for PoolInputTest
process TEST = {
	path p = {Analysis}
	module Analysis = OtherThingAnalyzer {untracked int32 debugLevel = 1}
	source = PoolRASource {
		untracked vstring fileNames = {
			'file:${LOCAL_TMP_DIR}/PoolInputTest.root',
			'PoolTest.root',
			'file:${LOCAL_TMP_DIR}/PoolInputOther.root'
		}
		untracked string catalog = '${LOCAL_TMP_DIR}/PoolInputTestCatalog.xml'
		untracked int32 maxEvents = -1
	}
}
!
cmsRun --parameter-set ${LOCAL_TMP_DIR}/PoolInputTest.cfg || die 'Failure using PoolInputTest.cfg' $?

