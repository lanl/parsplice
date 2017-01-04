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


#include "ldds.h"
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>

BDBLocalDataStore::BDBLocalDataStore(){

    
};

BDBLocalDataStore::~BDBLocalDataStore(){
    
    //close the databases
    for(std::map<unsigned int, Db*>::iterator it=dbm.begin();it!=dbm.end();it++){
        it->second->close(0);
        delete it->second;
    }
    
    delete env;
};

int BDBLocalDataStore::initialize(std::string homeDir_, std::string baseName_, unsigned long scratchSize_){
    
    BOOST_LOG_SEV(lg,info) <<"BDBLocalDataStore::initialize "<< homeDir_<<" "<<baseName_<<" "<<scratchSize_<<std::endl;

    
    homeDir=homeDir_;
    baseName=baseName_;
    scratchSize=scratchSize_;
    
    boost::filesystem::path p(homeDir);
    boost::filesystem::create_directories( p.parent_path().string() );

    
    //initialize the environment
    u_int32_t envFlags=
    DB_CREATE        |  // Create the environment if it does not exist
    //DB_PRIVATE |
    DB_RECOVER     |  // Run normal recovery.
    //DB_INIT_LOCK   |  // Initialize the locking subsystem
    DB_INIT_LOG    |  // Initialize the logging subsystem
    DB_INIT_TXN    |  // Initialize the transactional subsystem. This
    DB_AUTO_COMMIT |
    DB_INIT_MPOOL    |  // Initialize the memory pool (in-memory cache)
    DB_THREAD;          // Cause the environment to be free-threaded
    

    try {
        // create and open the environment
        env = new DbEnv(0);
        env->set_lk_detect(DB_LOCK_MINWRITE);
        
        env->set_error_stream(&std::cerr);
        env->set_cachesize(scratchSize, 0, 0);
        env->open(homeDir.c_str(), envFlags, 0644);
    }
    catch(DbException &e) {
        std::cerr << "Error opening database environment: "
        << homeDir << std::endl;
        std::cerr << e.what() << std::endl;
        return (EXIT_FAILURE);
    }
    return EXIT_SUCCESS;
};


int BDBLocalDataStore::put(unsigned int dbKey, Rd &key, Rd &data){
    
    
    //BOOST_LOG_SEV(lg,info) <<"#BDBLocalDataStore::put  count "<<count(dbKey,key)<<" "<<dbType[dbKey];
    BOOST_LOG_SEV(lg,debug) <<"#BDBLocalDataStore::put ";
    logkd(lg,dbKey,key,data);
    
    
    if(dbType.count(dbKey)==0){
        return DBKEY_NOTFOUND;
    }
    //convert raw data into Dbt's
    Dbt k(&(key[0]), u_int32_t(key.size()*sizeof(Rdt)));
    Dbt d(&(data[0]), u_int32_t(data.size()*sizeof(Rdt)));
    
    u_int32_t flags = 0 ;
    if(dbType[dbKey]){
        flags =  DB_NOOVERWRITE;
    }
    
    ml.lock();
    //put
    u_int32_t status=dbm[dbKey]->put(NULL, &k, &d, flags);
    ml.unlock();
    
    BOOST_LOG_SEV(lg,trace) <<"#BDBLocalDataStore::put  count "<<count(dbKey,key);
    BOOST_LOG_SEV(lg,trace) <<"#BDBLocalDataStore::put status "<<status;
    
    return status;
};

int BDBLocalDataStore::get(unsigned int dbKey, Rd &key, Rd &data){
    //BOOST_LOG_SEV(lg,info) <<"#BDBLocalDataStore::get count "<<count(dbKey,key);
    BOOST_LOG_SEV(lg,debug) <<"#BDBLocalDataStore::get ";
    
    data.clear();
    if( dbType.count(dbKey)==0 ){
        return DBKEY_NOTFOUND;
    }
    
    u_int32_t status=0;
    
    //convert raw data into Dbt's
    Dbt k(&(key[0]), u_int32_t(key.size()*sizeof(Rdt)));
    Dbt d;
    Dbc *cursorp;
    ml.lock();
    dbm[dbKey]->cursor(NULL, &cursorp, 0);
    status=cursorp->get(&k,&d,DB_SET);
    
    //std::cout<<status<<" "<<u_int32_t(DB_NOTFOUND)<<std::endl;
    if(status!=DB_NOTFOUND and status !=DB_KEYEMPTY){
        //convert Dbt into raw data
        //key.clear();
        key.resize( k.get_size()/sizeof(Rdt),0 );
        memcpy(&(key[0]),k.get_data(), k.get_size());
        
        data.clear();
        data.resize( d.get_size()/sizeof(Rdt),0 );
        memcpy(&(data[0]),d.get_data(), d.get_size());
        
        //delete the element we just read if it was not a named constant
        
        
        if(!dbType[dbKey]){
            cursorp->del(0);
        }
         
         
    }
    cursorp->close();
    ml.unlock();
    logkd(lg,dbKey,key,data);
    //BOOST_LOG_SEV(lg,info) <<"#BDBLocalDataStore::get count "<<count(dbKey,key);
    //BOOST_LOG_SEV(lg,info) <<"#BDBLocalDataStore::get status "<<status;
    
    return status;
};


unsigned int BDBLocalDataStore::count(unsigned int dbKey, Rd &key){
    
    if( dbType.count(dbKey)==0 ){
        return 0;
    }
    
    //convert raw data into Dbt's
    Dbt k(&(key[0]), u_int32_t(key.size()*sizeof(Rdt)));
    Dbt d;
    Dbc *cursorp;
    
    ml.lock();
    dbm[dbKey]->cursor(NULL, &cursorp, 0);
    u_int32_t status=cursorp->get(&k,&d,DB_SET);
    
    //BOOST_LOG_SEV(lg,info) <<"#BDBLocalDataStore::count status "<<status;
    
    //std::cout<<"STATUS: "<<status<<" "<<DB_KEYEMPTY<<" "<<DB_NOTFOUND<<std::endl;
    //if(status==DB_KEYEMPTY or status==DB_NOTFOUND){
    //    return 0;
    //}
    
    db_recno_t recno;
    u_int32_t status2;
    if(!(status==DB_KEYEMPTY or status==DB_NOTFOUND)){
        status2 =cursorp->count(&recno, 0);
    }
    
    BOOST_LOG_SEV(lg,trace) <<"#BDBLocalDataStore::count status "<<status2;
    cursorp->close();
    
    ml.unlock();
    
    if(status==DB_KEYEMPTY or status==DB_NOTFOUND){
        return 0;
    }
    return (unsigned int)(recno);
};




int BDBLocalDataStore::createDatabase(unsigned int dbKey,bool namedConstant){

    BOOST_LOG_SEV(lg,info) <<"BDBLocalDataStore::createDatabase "<< dbKey<<" "<<namedConstant<<" ";

    
    dbType[dbKey]=namedConstant;
    if( dbm.count(dbKey)==0 ){
        dbm[dbKey]=new Db(env,0);
        
        //create a database name
        std::string name=baseName+static_cast< std::ostringstream& >(std::ostringstream() << dbKey).str()+".db";
        //std::string name=baseName+(std::ostringstream() << dbKey).str()+".db";
        
        u_int32_t openFlags = DB_CREATE;  // Allow database creation
        dbm[dbKey]->set_flags(DB_DUP); //Allow duplicates
        
        //open the database
        int ret=dbm[dbKey]->open(NULL, name.c_str(), NULL, DB_HASH, openFlags, 0);
        BOOST_LOG_SEV(lg,debug) <<"#BDBLocalDataStore::createDatabase open "<<ret;
        
        
        /*
        DB_HASH_STAT *statp;
        ret=dbm[dbKey]->stat(NULL,statp, 0);
        
        
        BOOST_LOG_SEV(lg,debug) <<"#BDBLocalDataStore::createDatabase stat "<<ret;
        
        //BOOST_LOG_SEV(lg,debug) <<"#BDBLocalDataStore::createDatabase "<< (u_long) statp->hash_ndata<< "key/value pairs in the database";
        
        free(statp);
        */
    }

    return 0;
};

int BDBLocalDataStore::sync(){
    ml.lock();
    for(auto it=dbm.begin();it!=dbm.end();it++){
        it->second->sync(0);
    }
    ml.unlock();
    return 1;
};


STLLocalDataStore::STLLocalDataStore(){
    currentSize=0;
};

int STLLocalDataStore::initialize(std::string homeDir_, std::string baseName_, unsigned long maxSize_){
    maxSize=maxSize_*1000000000;
    return 0;
};




int STLLocalDataStore::put(unsigned int dbKey, Rd &key, Rd &data){
    
    BOOST_LOG_SEV(lg,debug) <<"#STLLocalDataStore::put ";
    
    
    if( dbType.count(dbKey)==0 ){
        BOOST_LOG_SEV(lg,error) <<"#STLLocalDataStore::put DBKEY_NOTFOUND";
        return DBKEY_NOTFOUND;
    }
    
    
    
    if( !dbType[dbKey] or (dbType[dbKey] and count(dbKey, key)==0) ){
        ml.lock();
        
        dbm[dbKey].insert(std::make_pair(key,data));
        currentSize+=data.size()*sizeof(Rdt);
        
        logkd(lg,dbKey,key,data);
        
        BOOST_LOG_SEV(lg,trace) <<"#STLLocalDataStore::put Inserted. Current size: "<<currentSize<<" of "<<maxSize;
        mru.touch( std::make_pair(dbKey,key ) );
        
        //we reached our maximum size. we need to prune
        while(currentSize>maxSize){
            
            BOOST_LOG_SEV(lg,trace) <<"#STLLocalDataStore::put Purging ";
    
            
            std::pair<bool, std::pair<unsigned int,Rd> > r=mru.back();
            if(r.first){
                //erase
                auto p=dbm[r.second.first].find(r.second.second);
                
                if(p!=dbm[r.second.first].end()){
                    currentSize-=p->second.size()*sizeof(Rdt);
                    dbm[r.second.first].erase(p);
                }
                if(dbm[r.second.first].count(r.second.second)==0){
                    mru.pop_back();
                }
                
                
                BOOST_LOG_SEV(lg,trace) <<"#STLLocalDataStore::put Purging. Current size: "<<currentSize<<" of "<<maxSize;

            }
            else{
                break;
            }
        }
        ml.unlock();
    }
    BOOST_LOG_SEV(lg,debug) <<"#STLLocalDataStore::put done";
    
    return 0;

};

int STLLocalDataStore::get(unsigned int dbKey, Rd &key, Rd &data){
    BOOST_LOG_SEV(lg,debug) <<"#STLLocalDataStore::get";
    data.clear();
   
    if( dbType.count(dbKey)==0 ){
        BOOST_LOG_SEV(lg,error) <<"#STLLocalDataStore::get DBKEY_NOTFOUND";
        return DBKEY_NOTFOUND;
    }
    
    ml.lock();
    if(dbm[dbKey].count(key)==0){
        ml.unlock();
        BOOST_LOG_SEV(lg,error) <<"#STLLocalDataStore::get DB_NOTFOUND";
        return DB_NOTFOUND;
    }
    

    data=dbm[dbKey].find(key)->second;
    logkd(lg,dbKey,key,data);
        
    mru.touch(std::make_pair(dbKey,key));
    
    if(!dbType[dbKey]){
        //need to erase from the db
        auto p=dbm[dbKey].find(key);
        currentSize-=p->second.size()*sizeof(Rdt);
        dbm[dbKey].erase(p);
        mru.erase(std::make_pair(dbKey,key));
    }
    ml.unlock();
    BOOST_LOG_SEV(lg,debug) <<"#STLLocalDataStore::get done";
    BOOST_LOG_SEV(lg,trace) <<"#STLLocalDataStore::get Current size: "<<currentSize<<" of "<<maxSize;
   
    return 0;

};

unsigned int STLLocalDataStore::count(unsigned int dbKey, Rd &key){
    ml.lock();
    unsigned int c;
    if( dbType.count(dbKey)==0 ){
        c=0;
    }
    else{
        c=dbm[dbKey].count(key);
    }
    ml.unlock();
    return c;
};

int STLLocalDataStore::createDatabase(unsigned int dbKey, bool namedConstant){
    
    if(dbm.count(dbKey)==0){
        dbm[dbKey]=std::multimap< Rd, Rd>();
        dbType[dbKey]=namedConstant;
        dbSize[dbKey]=0;
    }
    
    return 0;
};

int STLLocalDataStore::sync(){
    return 1;
};

