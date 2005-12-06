#include "DQMServices/WebComponents/interface/ContentReader.h"


// fills the list with all subdirectories of dir containing something
void ContentReader::give_subdirs(std::string dir, std::list<std::string> &subdirs, std::string mode)
{

  std::string old = mui->pwd();
  mui->setCurrentFolder(dir);

  std::vector<std::string> all_subdirs = mui->getSubdirs();

  std::cout << "The directory " << mui->pwd() << " contains " << all_subdirs.size() << " subdirectories!" << endl;

  subdirs.push_back("..");
  for (unsigned int i = 0; i < all_subdirs.size(); i++)
    {
      std::string subdir = dir + all_subdirs[i];
      if (mode == "User")
	{
	  if (mui->containsAnyMonitorable(subdir))
	    {
	      subdirs.push_back(all_subdirs[i]);
	    }
	}
      else if (mode == "SuperUser")
	{
	  subdirs.push_back(all_subdirs[i]);
	}
    }

  mui->setCurrentFolder(old);
}

// fills the list with existing files in dir
void ContentReader::give_files(std::string dir, std::list<std::string> &files, bool only_contents)
{
  cout << "ContentReader was asked to find the subfiles of : " << dir << endl;

  // for the interesting directory
  std::string old = mui->pwd();
  mui->setCurrentFolder(dir);
  // get the monitorable files
  std::vector<std::string> all_files = mui->getMEs();  

  if (only_contents)
    {
      // and copy those that exist in the files list
      for (unsigned int i = 0; i < all_files.size(); i++)
	{
	  std::string filename = mui->pwd() + "/" + all_files[i];
	  // if it does exist
	  if (mui->get(filename) != NULL) 
	    {
	      files.push_back(filename);
	    }
	}
      // go back to original folder:  
      mui->setCurrentFolder(old);
    }
  else
    {
      // and copy those that exist in the files list
      for (unsigned int i = 0; i < all_files.size(); i++)
	{
	  std::string filename = mui->pwd() + "/" + all_files[i];
	  files.push_back(filename);
	}
      // go back to original folder:  
      mui->setCurrentFolder(old);

    }
}


// returns a pointer to a requested filename, or NULL if not available
MonitorElement * ContentReader::give_ME(std::string filename)
{
  MonitorElement *requested = mui->get(filename);
  return requested;
}
