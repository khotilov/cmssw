#ifndef SISTRIPCORALIFACE_H
#define SISTRIPCORALIFACE_H


#include "CondCore/DBCommon/interface/RelationalStorageManager.h"
#include "CondCore/DBCommon/interface/AuthenticationMethod.h"
#include "CondCore/DBCommon/interface/SessionConfiguration.h"
#include "CondCore/DBCommon/interface/ConnectionConfiguration.h"
#include "CondCore/DBCommon/interface/MessageLevel.h"
#include "CondCore/DBCommon/interface/DBSession.h"
#include "CondCore/DBCommon/interface/Exception.h"
#include "CoralBase/TimeStamp.h"

#include <iterator>
#include <iostream>
#include <string>
#include <map>


	class SiStripCoralIface
	{
		public:	
			SiStripCoralIface(std::string connectionString, std::string authenticationPath);
			virtual ~SiStripCoralIface();
			void doQuery(coral::TimeStamp startTime, coral::TimeStamp endTime, std::vector<coral::TimeStamp>&, std::vector<uint32_t>&,std::vector<uint32_t>&  );
		private:
			void initialize();

			std::string m_connect;

			std::map<std::string,unsigned int> m_id_map;
			cond::DBSession* session;
			cond::RelationalStorageManager* m_coraldb;
	};
#endif
