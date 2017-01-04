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


#ifndef npdds_hpp
#define npdds_hpp

#include <stdio.h>


#include "dds.h"
#include "ParSpliceCommon.h"
#include "ParSpliceTypes.h"

#include <chrono>
#include <map>
#include <set>
#include <vector>
#include <list>
#include <queue>
#include <pthread.h>
#include <unistd.h>
#include <memory>
#include <atomic>
#include <future>

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>

#define DDS_PUT 1 << 1
#define DDS_GET 1 << 2
#define DDS_FORWARD 1 << 3



#define NBUFFERSERV 16


int mpiPack( Rdt *buffer, unsigned int &dbKey, Rd &key, Rd &data, MPI_Comm comm);
void mpiUnpack( Rdt *buffer, unsigned int &dbKey, Rd &key, Rd &data, MPI_Comm comm);


//abstract distributed data-store class

class AbstractDDS{
public:
    virtual void put(int mode, int dbKey, Rd &key, Rd &data)=0;
    virtual std::shared_future<Rd> get(int mode, int dbKey, Rd &key)=0;
    virtual void put(int mode, int dbKey, unsigned key, Rd &data)=0;
    virtual std::shared_future<Rd> get(int mode, int dbKey, unsigned key)=0;
    
    virtual bool serverSingle()=0;
    virtual int sync()=0;
};


template <class T> class dds : public AbstractDDS {
public:
    
    dds(MPI_Comm comm_, int parent_, std::string homeDir_, std::string baseName_,unsigned long maxSize_, std::map<unsigned int, bool>& dbType_, bool threadMultiple_) {
        BOOST_LOG_SEV(lg,debug) <<"#DDS::dds ";
        
        comm=comm_;
        parent=parent_;
        maxSize=maxSize_;
        homeDir=homeDir_;
        baseName=baseName_;
        dbType=dbType_;
        
        
        //create  directories
        boost::filesystem::path p(homeDir);
        boost::filesystem::create_directories( p.parent_path().string() );
        
        lstore.initialize(homeDir,baseName, maxSize);
        for(auto it=dbType.begin();it!=dbType.end();it++){
            lstore.createDatabase(it->first, it->second);
        }
        
        commBufferPut=std::vector<Rdt>(BUFFER_SIZE,0);
        commBufferGet=std::vector<Rdt>(BUFFER_SIZE,0);
        //commBufferServ=std::vector<Rdt>(BUFFER_SIZE,0);
        
        //multi-buffer for non-blocking receives
        commBufferServv= std::vector< std::vector<Rdt> >(NBUFFERSERV);
        for(int i=0;i<NBUFFERSERV;i++){
            commBufferServv[i]=std::vector<Rdt>(BUFFER_SIZE,0);
        }
        
        threadMultiple=threadMultiple_;
        if(threadMultiple){
            serverThread=std::thread(&dds::server, this);
        }
        else{
            completed=0;
            //post an initial non-blocking receive
            //MPI_Irecv(&(commBufferServ[0]),BUFFER_SIZE,MPI_PACKED,MPI_ANY_SOURCE,MPI_ANY_TAG,comm,&r);
            for(int i=0;i<NBUFFERSERV;i++){
                MPI_Irecv(&(commBufferServv[i][0]),BUFFER_SIZE,MPI_PACKED,MPI_ANY_SOURCE,MPI_ANY_TAG,comm,&rv[i]);
            }
        }

    };
    
    

    
    
    ~dds(){
        serverThread.join();
    };
    
    
    bool process( Request &req){
        
        bool success;
        
        BOOST_LOG_SEV(lg,debug) <<"#DDS::process";
        
        bool put = req.type & DDS_PUT;
        bool get= req.type & DDS_GET;
        
        
        if(put){
            success=processPut(req);
        }
        if(get){
            success=processGet(req);
        }
        
        return success;
    };
    
    bool processPut(Request &req){
        //std::vector<Rdt> commBuffer(BUFFER_SIZE,0);
        unsigned *p=reinterpret_cast<unsigned*>(&(req.key[0]));
        BOOST_LOG_SEV(lg,info) <<"#DDS::processPut "<<req.type<<" "<<*p;
        
        //if the requested database does not exist, ignore the request
        if(dbType.count(req.dbKey)==0){
            BOOST_LOG_SEV(lg,error) <<"#DDS::processPut Database " <<req.dbKey<<" does not exist";
            return false;
        }
        
        //is the requested data already in the local store?
        bool inLocalStore=lstore.count(req.dbKey,req.key);
        //does this database store named constants?
        bool namedConstant=dbType[req.dbKey];
        //should we add duplicate entries?
        bool duplicate=!namedConstant;
        
        BOOST_LOG_SEV(lg,trace) <<"#DDS::processPut inLocalStore " <<lstore.count(req.dbKey,req.key)<<" namedConstant "<<namedConstant<<" duplicate "<<duplicate;
        
        
        //should we forward the request?
        bool forward= req.type & DDS_FORWARD;
        //BOOST_LOG_SEV(lg,info) <<"#DDS::processPut forward "<<forward;
        //only forward if we are not at the root
        forward = forward and (parent>=0);
        //BOOST_LOG_SEV(lg,info) <<"#DDS::processPut forward "<<forward;
        //dont forward if this is a named constant that we already forwarded in the past
        forward = forward and ( !namedConstant or (namedConstant and forwardedPuts.count(std::make_pair(req.dbKey,req.key))==0)  );
        //BOOST_LOG_SEV(lg,info) <<"#DDS::processPut forward "<<forward;
        //don't forward if the request came from that parent
        forward = forward and !(req.source==parent);
        //BOOST_LOG_SEV(lg,info) <<"#DDS::processPut forward "<<forward;
        
        //should we even store locally?
        bool store=!(forward and not namedConstant);
        
        
        if(forward){
            BOOST_LOG_SEV(lg,debug) <<"#DDS::processPut Forwarding to parent "<<parent;
            int size=mpiPack( &(commBufferPut[0]), req.dbKey, req.key, req.data, comm );
            
            //start a non-blocking send
            MPI_Request r;
            MPI_Status s;
            MPI_Issend(&(commBufferPut[0]),size,MPI_PACKED, parent, req.type, comm , &r);
            
            std::chrono::high_resolution_clock::time_point start=std::chrono::high_resolution_clock::now();
            std::chrono::milliseconds timeout=std::chrono::milliseconds(100);
            int completed=0;
            int cancelled=0;
            
            while(true){
                //the communication is taking too long. Attempt to cancel
                cancelled=0;
                completed=0;
                if( std::chrono::high_resolution_clock::now() - start > timeout ){
                     BOOST_LOG_SEV(lg,error) <<"#DDS::processPut initiating cancel of forward request  "<<parent;
                    //cancel
                    MPI_Cancel(&r);
                    //wait for cancellation or completion
                    MPI_Wait(&r, &s);
                    //did the cancellation succeed?
                    MPI_Test_cancelled(&s, &cancelled);
                    if(!cancelled){
                        completed=1;
                    }
                    if(cancelled){
                        completed=0;
                    }
                    BOOST_LOG_SEV(lg,error) <<"#DDS::processPut cancelling forward request  "<<parent;
                }
                else{
                    //test whether the communication was successful
                    MPI_Test(&r, &completed, &s);
                }
                
                if(completed){
                    forwardedPuts.insert(std::make_pair(req.dbKey,req.key));
                    BOOST_LOG_SEV(lg,debug) <<"#DDS::processPut request forwarded "<<parent;
                    break;
                }
                if(cancelled){
                    BOOST_LOG_SEV(lg,error) <<"#DDS::processPut forward request failed "<<parent;
                    return false;
                }
            }

            //MPI_Send(&(commBufferPut[0]),size,MPI_PACKED, parent, req.type, comm );
            //forwardedPuts.insert(std::make_pair(req.dbKey,req.key));
            //BOOST_LOG_SEV(lg,debug) <<"#DDS::processPut request forwarded "<<parent;
        }
        
        if(store){
            //add to local store
            if(duplicate or !inLocalStore){
                BOOST_LOG_SEV(lg,trace) <<"#DDS::processPut Adding to local store";
                logRequest(lg,req);
                lstore.put(req.dbKey,req.key, req.data);
            }
        }
        
        //if we kept a local copy, a put can allow us to fulfill pending get requests
        BOOST_LOG_SEV(lg,debug) <<"#DDS::processPut Processing follow-on requests "<<pendingRequests.count(std::make_pair(req.dbKey,req.key) );
        
        std::list< Request > requestsToProcess;
        
        if(store){
            if(namedConstant){
                //we can process all pending requests for this item
                auto p = pendingRequests.equal_range( std::make_pair(req.dbKey,req.key) );
                for(auto it=p.first;it!=p.second;it++){
                    requestsToProcess.push_back(it->second);
                }
                pendingRequests.erase( std::make_pair(req.dbKey,req.key) );
            }
            else{
                //we can process only one pending request for this item
                if(pendingRequests.count( std::make_pair(req.dbKey,req.key) ) >0){
                    auto it=pendingRequests.find( std::make_pair(req.dbKey,req.key) );
                    requestsToProcess.push_back(it->second);
                    pendingRequests.erase( it );
                }
            }
            
            //process pending requests
            for(auto it=requestsToProcess.begin(); it!=requestsToProcess.end(); it++){
                BOOST_LOG_SEV(lg,trace) <<"#DDS::processPut Process request "<<it->type<<" from "<<it->source;
                processGet(*it);
            }
        }
        BOOST_LOG_SEV(lg,debug) <<"#DDS::processPut Done";
        //pthread_mutex_unlock(&processLock);
        
        return true;
    };
    
    bool processGet(Request &req){
        //std::vector<Rdt> commBuffer(BUFFER_SIZE,0);
        //pthread_mutex_lock(&processLock);
        unsigned *p=reinterpret_cast<unsigned*>(&(req.key[0]));
        BOOST_LOG_SEV(lg,info) <<"#DDS::processGet "<<req.type<<" "<<*p;
        
        //this is a get request, so data shoud be empty anyway
        //clear it up to be sure
        req.data.clear();
        
        //if the requested database does not exist, ignore
        if(dbType.count(req.dbKey)==0 ){
            BOOST_LOG_SEV(lg,error) <<"#DDS::processGet Database "<<req.dbKey<<" does not exist";
            return false;
        }
        
        //is the requested data in the local store?
        bool inLocalStore=lstore.count(req.dbKey,req.key);
        
        //does this database store named constants?
        bool namedConstant=dbType[req.dbKey];
        
        BOOST_LOG_SEV(lg,debug) <<"#DDS::processGet inLocalStore: "<<lstore.count(req.dbKey,req.key)<<" namedConstant "<<namedConstant;
        
        //Answer the request
        if(inLocalStore){
            
            int s=lstore.get(req.dbKey,req.key,req.data);
            logRequest(lg,req);
            
            BOOST_LOG_SEV(lg,debug) <<"#DDS::processGet sending reply to "<<req.source;
            //this was an mpi request
            if( req.source>=0){
                int size=mpiPack( &(commBufferGet[0]), req.dbKey, req.key, req.data,comm);
                
                
                
                //start a non-blocking send
                MPI_Request r;
                MPI_Status s;
                MPI_Issend(&(commBufferGet[0]),size,MPI_PACKED, req.source, DDS_PUT | DDS_FORWARD,comm, &r);
                
                std::chrono::high_resolution_clock::time_point start=std::chrono::high_resolution_clock::now();
                std::chrono::milliseconds timeout=std::chrono::milliseconds(100);
                int completed=0;
                int cancelled=0;
                while(true){
                    completed=0;
                    cancelled=0;
                    if( std::chrono::high_resolution_clock::now() - start > timeout ){
                        MPI_Cancel(&r);
                        MPI_Wait(&r, &s);
                        MPI_Test_cancelled(&s, &cancelled);
                        if(!cancelled){
                            completed=1;
                        }
                        if(cancelled){
                            completed=0;
                        }
                        BOOST_LOG_SEV(lg,error) <<"#DDS::processGet concelling reply   "<<req.source;
                    }
                    else{
                        MPI_Test(&r, &completed, &s);
                    }
                    
                    
                    if(completed){
                        BOOST_LOG_SEV(lg,debug) <<"#DDS::processGet reply sent  "<<req.source;
                        break;
                    }
                    if(cancelled){
                        BOOST_LOG_SEV(lg,error) <<"#DDS::processGet reply failed "<<req.source;
                        return false;
                    }
                    
                }
                
                
                
                //MPI_Send(&(commBufferGet[0]),size,MPI_PACKED, req.source, DDS_PUT | DDS_FORWARD,comm);
                //BOOST_LOG_SEV(lg,debug) <<"#DDS::processGet reply sent "<<req.source;
            }
            //this was a local request
            else{
                //find the corresponding pending promise
                promiseLock.lock();
                if(pendingPromises.count(std::make_pair(req.dbKey, req.key) )> 0){
                    auto it=pendingPromises.find(std::make_pair(req.dbKey, req.key));
                    it->second.set_value(req.data);
                    //delete the promise
                    pendingPromises.erase(it);
                }
                promiseLock.unlock();
            }
        }
        else{
            BOOST_LOG_SEV(lg,trace) <<"#DDS::processGet saving pending request "<<req.source;
            //if not, keep the request in store to fulfill later
            pendingRequests.insert( std::make_pair(std::make_pair(req.dbKey,req.key),req) );
        }
        
        bool forward= req.type & DDS_FORWARD;
        //BOOST_LOG_SEV(lg,info) <<"#DDS::processGet forward "<<forward;
        //dont forward if this is a named constant we already requested
        //forward = forward and ( !namedConstant or (namedConstant and pendingRequests.count(std::make_pair(req.dbKey,req.key))==0)  );
        forward = forward and ( !namedConstant or (namedConstant and forwardedGets.count(std::make_pair(req.dbKey,req.key))==0)  );
        //BOOST_LOG_SEV(lg,info) <<"#DDS::processGet forward "<<forward;
        //only forward if we are not a root
        forward = forward and (parent>=0);
        //BOOST_LOG_SEV(lg,info) <<"#DDS::processGet forward "<<forward;
        //dont forward if we fulfilled the request
        forward = forward and !inLocalStore;
        //BOOST_LOG_SEV(lg,info) <<"#DDS::processGet forward "<<forward;
        
        if(forward){
            BOOST_LOG_SEV(lg,debug) <<"#DDS::processGet Forwarding to parent "<<parent;
            unsigned int size=mpiPack( &(commBufferGet[0]), req.dbKey, req.key, req.data, comm);
            
            //start a non-blocking send
            MPI_Request r;
            MPI_Status s;
            MPI_Issend(&(commBufferGet[0]),size,MPI_PACKED, parent, req.type,comm, &r);
            
            std::chrono::high_resolution_clock::time_point start=std::chrono::high_resolution_clock::now();
            std::chrono::milliseconds timeout=std::chrono::milliseconds(100);
            int completed=0;
            int cancelled=0;
            
            
            while(true){
                completed=0;
                cancelled=0;
                if( std::chrono::high_resolution_clock::now() - start > timeout ){
                    MPI_Cancel(&r);
                    MPI_Wait(&r, &s);
                    MPI_Test_cancelled(&s, &cancelled);
                    if(!cancelled){
                        completed=1;
                    }
                    if(cancelled){
                        completed=0;
                    }
                    BOOST_LOG_SEV(lg,error) <<"#DDS::processGet request cancelled   "<<parent;
                }
                else{
                    MPI_Test(&r, &completed, &s);
                }
                
                if(completed){
                    BOOST_LOG_SEV(lg,debug) <<"#DDS::processGet request forwarded  "<<parent;
                    break;
                }
                if(cancelled){
                    BOOST_LOG_SEV(lg,error) <<"#DDS::processGet forward failed "<<parent;
                    return false;
                }
                
            }
            
            if(completed){
                BOOST_LOG_SEV(lg,trace) <<"#DDS::processGet saving pending request "<<req.source;
                forwardedGets.insert(std::make_pair(req.dbKey,req.key));
            }
            //MPI_Send(&(commBufferGet[0]),size,MPI_PACKED, parent, req.type,comm);
            //BOOST_LOG_SEV(lg,debug) <<"#DDS::processGet request forwarded "<<parent;
        }
        else{
            
        }
        BOOST_LOG_SEV(lg,debug) <<"#DDS::processGet Done";
        return true;
    };
    
    static void *serverHandle(void *context){
        dds<T> *p=static_cast< dds<T> * >(context);
        p->server();
        return NULL;
    };
    
    void server(){
        //std::vector<Rdt> commBuffer(BUFFER_SIZE,0);
        
                
        BOOST_LOG_SEV(lg,debug) <<"#DDS::server ";
        
        completed=0;
        //post an initial non-blocking receive
        //MPI_Irecv(&(commBufferServ[0]),BUFFER_SIZE,MPI_PACKED,MPI_ANY_SOURCE,MPI_ANY_TAG,comm,&r);
        for(int i=0;i<NBUFFERSERV;i++){
            MPI_Irecv(&(commBufferServv[i][0]),BUFFER_SIZE,MPI_PACKED,MPI_ANY_SOURCE,MPI_ANY_TAG,comm,&rv[i]);
        }
        while(true){
            for(int i=0;i<NBUFFERSERV;i++){
                completed=0;
                MPI_Test(&rv[i],&completed,&status);
                
                //std::cout<<completed<<std::endl;
                
                if(completed){
                    Rd key;
                    Rd data;
                    unsigned int dbKey;
                    mpiUnpack( &(commBufferServv[i][0]), dbKey, key, data,comm);
                    
                    int task=status.MPI_TAG;
                    int source=status.MPI_SOURCE;
                    
                    BOOST_LOG_SEV(lg,trace) <<"#DDS::server Received request "<<task<<" from "<<source;
                    
                    
                    Request req(task, dbKey,key, data,source);
                    logRequest(lg,req);
                    
                    requestLock.lock();
                    incomingRequests.push(req);
                    requestLock.unlock();
                    
                    //post the next receive
                    completed=0;
                    MPI_Irecv(&(commBufferServv[i][0]),BUFFER_SIZE,MPI_PACKED,MPI_ANY_SOURCE,MPI_ANY_TAG,comm,&rv[i]);
                    
                }
            }
            //process pending requests
            requestLock.lock();
            std::queue<Request> failedRequests;
            while(incomingRequests.size()>0){
                logRequest(lg,incomingRequests.front());
                bool success=process(incomingRequests.front());
                if(!success){
                    failedRequests.push(incomingRequests.front());
                }
                incomingRequests.pop();
            }
            while(failedRequests.size()>0){
                incomingRequests.push(failedRequests.front());
                failedRequests.pop();
            }
            requestLock.unlock();
        }
    };
    
    
    bool serverSingle(){
        
        if(threadMultiple){
            //no need to do anything if there is already a server thread running
            return false;
        }
        
        //BOOST_LOG_SEV(lg,debug) <<"#DDS::serverSingle ";
        bool messageReceived=false;
        
        for(int i=0;i<NBUFFERSERV;i++){
            completed=0;
            MPI_Test(&rv[i],&completed,&status);
            
            if(completed){
                messageReceived=true;
                Rd key;
                Rd data;
                unsigned int dbKey;
                //mpiUnpack( &(commBufferServ[0]), dbKey, key, data,comm);
                mpiUnpack( &(commBufferServv[i][0]), dbKey, key, data,comm);
                
                int task=status.MPI_TAG;
                int source=status.MPI_SOURCE;
                
                BOOST_LOG_SEV(lg,debug) <<"#DDS::server Received request "<<task<<" from "<<source<<" in buffer "<<i;
                
                
                Request req(task, dbKey,key, data,source);
                logRequest(lg,req);
                
                requestLock.lock();
                incomingRequests.push(req);
                requestLock.unlock();
                
                //post the next receive
                completed=0;
                //MPI_Irecv(&(commBufferServ[0]),BUFFER_SIZE,MPI_PACKED,MPI_ANY_SOURCE,MPI_ANY_TAG,comm,&r);
                MPI_Irecv(&(commBufferServv[i][0]),BUFFER_SIZE,MPI_PACKED,MPI_ANY_SOURCE,MPI_ANY_TAG,comm,&rv[i]);
            }
            
            //process pending requests
            requestLock.lock();
            std::queue<Request> failedRequests;
            while(incomingRequests.size()>0){
                logRequest(lg,incomingRequests.front());
                bool success=process(incomingRequests.front());
                if(!success){
                    failedRequests.push(incomingRequests.front());
                }
                incomingRequests.pop();
            }
            while(failedRequests.size()>0){
                incomingRequests.push(failedRequests.front());
                failedRequests.pop();
            }
            requestLock.unlock();
        }
        return false;
        //return messageReceived;
    };

    
    
    
    
    virtual void put(int mode, int dbKey, unsigned key, Rd &data){
        BOOST_LOG_SEV(lg,debug) <<"#DDS::put "<<mode<<" "<<dbKey<<" "<<key;
        Rd k(sizeof(unsigned)/sizeof(Rdt),0);
        unsigned *p=reinterpret_cast<unsigned*>(&(k[0]));
        *p=key;
        put(mode, dbKey, k, data);
    };
    
    
    virtual std::shared_future<Rd> get(int mode, int dbKey, unsigned key){
        BOOST_LOG_SEV(lg,debug) <<"#DDS::get "<<mode<<" "<<dbKey<<" "<<key;
        Rd k(sizeof(unsigned)/sizeof(Rdt),0);
        unsigned *p=reinterpret_cast<unsigned*>(&(k[0]));
        *p=key;
        return get(mode, dbKey, k);
    };
    
    virtual void put(int mode, int dbKey, Rd &key, Rd &data){
        //BOOST_LOG_SEV(lg,debug) <<"#DDS::put "<<mode<<" "<<dbKey<<" "<<key;
        Request req(mode, dbKey,key, data, -1);
        logRequest(lg,req);
        
        
        requestLock.lock();
        incomingRequests.push(req);
        requestLock.unlock();
    };
    
    
    virtual std::shared_future<Rd> get(int mode, int dbKey, Rd &key){
        Rd data;
        Request req(mode, dbKey,key, data, -1);
        logRequest(lg,req);
        
        requestLock.lock();
        incomingRequests.push(req);
        requestLock.unlock();
        
        std::promise<Rd> p;
        promiseLock.lock();
        auto it=pendingPromises.insert(std::make_pair(std::make_pair(dbKey,key), std::promise<Rd>() ));
        std::shared_future<Rd> sf(it->second.get_future());
        promiseLock.unlock();
    
        return sf;
    };
    
    
    virtual int sync(){
        BOOST_LOG_SEV(lg,debug) <<"#DDS::sync ";
        return lstore.sync();
    };
    
private:
    
    //local store
    T lstore;
    
    std::vector<Rdt> commBufferPut;
    std::vector<Rdt> commBufferGet;
    //std::vector<Rdt> commBufferServ;
    
     std::vector< std::vector<Rdt> > commBufferServv;
    
    
    std::queue<Request> incomingRequests;
    
    //the keys of the named constants put that we already forwarded to our parent
    std::set< std::pair<unsigned int, Rd> > forwardedPuts;
    //the keys of the named constants get that we already requested to our parent
    std::set< std::pair<unsigned int, Rd> > forwardedGets;
    
    
    //a log of the pending (get) requests
    std::multimap< std::pair<unsigned int, Rd>, Request > pendingRequests;
    
    std::multimap< std::pair<unsigned int, Rd>, std::promise<Rd> > pendingPromises;
    
    std::map<unsigned int, bool> dbType;
    
    MPI_Comm comm;
    int parent;
    
    
    bool (*inferType)(unsigned int);
    std::string homeDir;
    std::string baseName;
    unsigned long maxSize;
    
    boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg;
    
    
    std::thread serverThread;
    parsplice::SpinLock promiseLock;
    parsplice::SpinLock requestLock;
    
    
    bool threadMultiple;
    
    //MPI_Request r;
    MPI_Request rv[NBUFFERSERV];
    MPI_Status status;
    int completed;
    
    
    
    std::chrono::milliseconds checkpointDelay;
    
};





#endif /* npdds_hpp */
