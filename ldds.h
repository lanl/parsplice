/*
 Copyright (c) 2016, Los Alamos National Security, LLC
 All rights reserved.
 Copyright 2016. Los Alamos National Security, LLC. This software was produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from LANL.
 
 Additionally, redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 1.      Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 2.      Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 3.      Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL, the U.S. Government, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


#ifndef dds_dds_h
#define dds_dds_h


#include <iostream>
#include <db_cxx.h>

#include <pthread.h>
#include <unistd.h>
#include <sstream>

#include <map>
#include <vector>
#include <string>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/filesystem.hpp>


#include "dds.h"
//#include "ParSpliceTypes.h"

using namespace boost::multi_index;

//template <typename Item>
class SpinLock {
    std::atomic_flag locked = ATOMIC_FLAG_INIT ;
public:
    void lock() {
        while (locked.test_and_set(std::memory_order_acquire)) { ; }
    }
    void unlock() {
        locked.clear(std::memory_order_release);
    }
};

class MRU
{
    typedef boost::multi_index::multi_index_container<
    std::pair<unsigned int, Rd> ,
    boost::multi_index::indexed_by<
    boost::multi_index::sequenced<>,
    boost::multi_index::hashed_unique<boost::multi_index::identity<std::pair<unsigned int, Rd> > >
    >
    > item_list;
    
   
public:
    //typedef Item                         item_type;
    typedef typename item_list::iterator iterator;
    std::size_t max_num_items;
    

    void touch( std::pair<unsigned int, Rd> item){
        std::pair<iterator,bool> p=il.push_front(item);
        
        if(!p.second){                     /* duplicate item */
            il.relocate(il.begin(),p.first); /* put in front */
        }
    }
    
    void erase( std::pair<unsigned int, Rd> item){
        il.get<1>().erase(item);
    }
    
    std::pair<bool, std::pair<unsigned int, Rd> > back(){
        std::pair<unsigned int, Rd> b;
        bool ok=il.size()>0;
        if(ok){
            b=il.back();
        }
        return std::make_pair(ok,b);
    }
    
    bool pop_back(){
        std::pair<unsigned int, Rd> b;
        bool ok=il.size()>0;
        if(ok){
            il.pop_back();
        }
        return ok;
    }
    
    
private:
    item_list   il;
    
};




//abstract base class for the local data store
class AbstractLocalDataStore{
public:
    virtual int put(unsigned int dbKey, Rd &key, Rd &data)=0;
    virtual int get(unsigned int dbKey, Rd &key, Rd &data)=0;
    virtual unsigned int count(unsigned int dbKey, Rd &key)=0;
    virtual int createDatabase(unsigned int dbKey, bool namedConstant)=0;
    virtual int initialize(std::string homeDir, std::string baseName, unsigned long maxSize)=0;
    virtual int sync()=0;
    
protected:
    boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg;
    
    SpinLock ml;
};

class BDBLocalDataStore : public AbstractLocalDataStore {
public:
    BDBLocalDataStore();
    ~BDBLocalDataStore();
    
    int initialize(std::string homeDir, std::string baseName, unsigned long maxSize);
    virtual int put(unsigned int dbKey, Rd &key, Rd &data);
    virtual int get(unsigned int dbKey, Rd &key, Rd &data);
    virtual unsigned int count(unsigned int dbKey, Rd &key);
    int createDatabase(unsigned int dbKey, bool namedConstant);
    int sync();
    
private:
    //int createDatabase(unsigned int dbKey);
    
    
    std::string homeDir;
    std::string baseName;
    unsigned long scratchSize;
    DbEnv *env;
    std::map<unsigned int, Db*> dbm;
    
    std::map<unsigned int, bool> dbType;
};



class STLLocalDataStore : public AbstractLocalDataStore {
public:
    STLLocalDataStore();
    //~STLLocalDataStore();
    
    int initialize(std::string homeDir, std::string baseName, unsigned long maxSize);
    virtual int put(unsigned int dbKey, Rd &key, Rd &data);
    virtual int get(unsigned int dbKey, Rd &key, Rd &data);
    virtual unsigned int count(unsigned int dbKey, Rd &key);
    int createDatabase(unsigned int dbKey, bool namedConstant);
    int sync();
private:
    unsigned int scratchSize;
    
    std::map<unsigned int, bool> dbType;
    std::map<unsigned int, unsigned int> dbSize;
    
    unsigned long maxSize;
    unsigned long currentSize;
    std::map< unsigned int, std::multimap< Rd, Rd> > dbm;
    
    MRU mru;
    
    
};

#endif
