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

#ifndef __ParSplice__NodeManager__
#define __ParSplice__NodeManager__

#include <stdio.h>
#include <set>
#include <chrono>

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>


#include "ParSpliceCommon.h"
#include "Validator.h"
#include "ParSpliceWorker.h"
#include "npdds.h"
#include "ldds.h"


#define NCOMMBUFF 16

class NodeManager{
    
public:
    
    /*
    NodeManager(MPI_Comm comm_, int rank_,std::set<int> children_,int parent_, int localBatchSize_, int globalBatchSize_, int nWorkers_, AbstractDDS &dds_) : dds(dds_){
        comm=comm_;
        children=children_;
        parent=parent_;
        nWorkers=nWorkers_;
        rank=rank_;
        
        nWork=nWorkers;
        nPrefetch=nWorkers;
        
        std::set<int> gclients=children;
        gclients.insert(rank);
        globalValidator.init(gclients,globalBatchSize_);
        
        std::set<int> lclients;
        for(int i=0;i<nWorkers;i++){
            lclients.insert(i);
        }
        localValidator.init(lclients,localBatchSize_);
        
        
        workers=std::vector<std::thread>(nWorkers);
        for(int i=0;i<nWorkers;i++){
            workers[i]=std::thread(&NodeManager::worker, this,i);
        }
        
        int rank=0;
        int nRanks=1;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &nRanks);
        
        std::string homeDir=str(boost::format("./r%i-hot/") % rank);
        std::string baseName="hot";
        unsigned long maxSize=1;
        
        hdds.initialize(homeDir, baseName, maxSize);
        hdds.createDatabase(1, false);
    };
     */
    
    
    NodeManager(MPI_Comm comm_, MPI_Comm dbComm_, int rank_,std::set<int> children_,int parent_, boost::property_tree::ptree &config, bool threadMultiple_){
        comm=comm_;
        children=children_;
        parent=parent_;
        rank=rank_;
        threadMultiple=threadMultiple_;
        
        start=std::chrono::high_resolution_clock::now();
        
        nWorkers=config.get<int>("ParSplice.Topology.WorkersPerNode");
        nWork=config.get<int>("ParSplice.NodeManager.WorkBufferLength");
        nPrefetch=config.get<int>("ParSplice.NodeManager.PrefetchBufferLength");
        
        checkpointDelay=std::chrono::milliseconds( config.get<unsigned>("ParSplice.Splicer.CheckpointDelay",100) );
        

        std::set<int> gclients=children;
        gclients.insert(rank);
        int globalBatchSize=config.get<int>("ParSplice.Validation.GlobalBatchSize",1);
        int validatorDelay=config.get<int>("ParSplice.Validation.ReleaseDelay",1);
        //globalValidator.init(gclients,validatorDelay);
        globalValidator.init(gclients,globalBatchSize);
        
        std::set<int> lclients;
        for(int i=0;i<nWorkers;i++){
            lclients.insert(i);
        }
        int localBatchSize=config.get<int>("ParSplice.Validation.LocalBatchSize",1);
        localValidator.init(lclients,localBatchSize);
        
        
        
        //create the databases
        std::map<unsigned int, bool> dbType;
        dbType[1]=true;
        
        std::string sdbType=config.get<std::string>("ParSplice.NodeManager.DB.StateDB.DBType");
        std::string base=config.get<std::string>("ParSplice.NodeManager.DB.StateDB.BaseName");
        std::string home=config.get<std::string>("ParSplice.NodeManager.DB.StateDB.Home");
        boost::trim(sdbType);
        boost::trim(base);
        boost::trim(home);
        unsigned long footprint=config.get<unsigned long>("ParSplice.NodeManager.DB.StateDB.Footprint");
        std::string rs=boost::str(boost::format("%1%" ) % rank );
        boost::replace_all(home, "__RANK__", rs);
        
        BOOST_LOG_SEV(lg,trace) <<"#NodeManager "<<parent<<" "<<home<<" "<<base<<" "<<footprint<<" ";
        
        if(sdbType.find("BDB")!=std::string::npos){
            sdb=new dds<BDBLocalDataStore>(dbComm_,parent,home,base,footprint,dbType,threadMultiple);
        }
        if(sdbType.find("STL")!=std::string::npos){
            sdb=new dds<STLLocalDataStore>(dbComm_,parent,home,base,footprint,dbType,threadMultiple);
        }
        
        
        std::string hdbType=config.get<std::string>("ParSplice.NodeManager.DB.HotDB.DBType");
        base=config.get<std::string>("ParSplice.NodeManager.DB.HotDB.BaseName");
        home=config.get<std::string>("ParSplice.NodeManager.DB.HotDB.Home");
        boost::trim(hdbType);
        boost::trim(base);
        boost::trim(home);
        footprint=config.get<unsigned long>("ParSplice.NodeManager.DB.HotDB.Footprint");
        boost::replace_all(home, "__RANK__", rs);
        BOOST_LOG_SEV(lg,trace) <<"#NodeManager "<<parent<<" "<<home<<" "<<base<<" "<<footprint<<" ";
        
        
        if(hdbType.find("BDB")!=std::string::npos){
            hdb=new BDBLocalDataStore;
        }
        if(hdbType.find("STL")!=std::string::npos){
            hdb=new STLLocalDataStore;
        }
        hdb->initialize(home, base, footprint);
        hdb->createDatabase(1, false);

        
        //hack to preserve processor affinity on SLURM systems
        drivers=NULL;
        drivers=new ParSpliceWorkerDriver[nWorkers];
        
        workers=std::vector<std::thread>(nWorkers);
        for(int i=0;i<nWorkers;i++){
            workers[i]=std::thread(&NodeManager::worker, this,i,std::ref(drivers[i]));
        }
    };

    
    
    
    ~NodeManager(){
        delete [] drivers;
        
    };
    
    
    void worker(int rank, ParSpliceWorkerDriver &driver){
        
        
        
        BOOST_LOG_SEV(lg,debug) <<"#NodeManager::worker";
        
        //initialize the driver
        //ParSpliceWorkerDriver driver;
        
        
        rand[rank].seed(rank*1234);
        
        
        while(true){
            BOOST_LOG_SEV(lg,trace) <<"#NodeManager::worker " <<rank<<" starting new segment";
            //get a random state from the work buffer
            workLock.lock();
            int nw=workBuffer.size();
            workLock.unlock();
            
            BOOST_LOG_SEV(lg,trace) <<"#NodeManager::worker number of elements in the work buffer "<<nw;
            
            if(nw>0){
                boost::random::uniform_int_distribution<int> d(0,nw-1);
                int iw=d(rand[rank]);
                workLock.lock();
                std::pair<parsplice::Label, Rd> w=workBuffer[iw];
                workLock.unlock();
                parsplice::Label wLabel=w.first;
                Rd wRef=w.second;
                
                BOOST_LOG_SEV(lg,info) <<"#NodeManager::worker "<<rank<<" running task "<<wLabel;
                
                //hot point
                Rd wHot;
                {
                    Rd k(sizeof(unsigned)/sizeof(Rdt),0);
                    unsigned *p=reinterpret_cast<unsigned*>(&(k[0]));
                    *p=wLabel;
                    hdb->get(1,k,wHot);
                    BOOST_LOG_SEV(lg,trace) <<"#NodeManager::worker hot point? "<<bool(wHot.size()>0);
                }
                parsplice::TrajDB seg;
                //run the driver
                //upon return, wLabel,wRef,wHot, are replaced by the new values
                
                parsplice::Trajectory t;
                parsplice::Label wLabel0=wLabel;
                
                std::chrono::high_resolution_clock::time_point segmentStart=std::chrono::high_resolution_clock::now();
                driver.generateSegment(wLabel, wRef, wHot, t);
                BOOST_LOG_SEV(lg,info) <<"#NodeManager::worker "<<rank<<" completed task "<<wLabel<<" in "<<std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now()-segmentStart).count()<< "s";
                
                seg[wLabel0].push_back(t);
                
                //BOOST_LOG_SEV(lg,info) <<"#NodeManager::worker "<<rank<<" completed task "<<wLabel;
                logTrajDB(lg,seg);
                
                //put the final state in the db
                sdb->put(DDS_PUT | DDS_FORWARD, 1, wLabel, wRef);
                //put the final hot point in the local store
                {
                    Rd k(sizeof(unsigned)/sizeof(Rdt),0);
                    unsigned *p=reinterpret_cast<unsigned*>(&(k[0]));
                    *p=wLabel;
                    hdb->put(1,k,wHot);
                }
                
                BOOST_LOG_SEV(lg,info) <<"#NodeManager::worker "<<rank<<" saved hot point "<<wLabel;
                
                //register the segments
                lSegmentLock.lock();
                localPendingSegments.push_back( std::make_pair<>(rank,seg) );
                lSegmentLock.unlock();
                
                BOOST_LOG_SEV(lg,info) <<"#NodeManager::worker "<<rank<<" logged task "<<wLabel;
            }
            else{
                BOOST_LOG_SEV(lg,debug) <<"#NodeManager::worker work buffer is empty";
                //std::this_thread::sleep_for(std::chrono::microseconds(1000));
            }
            
        }
        
        
       
        
        
        
        /*
        rand[rank].seed(rank*1234);
        
       
        boost::random::uniform_int_distribution<int> d(1,10);

        while(true){
            TrajDB db;
            Trajectory t;
            Visit v;
            
            
            //generate a segment
            v.label=d(rand[rank]);
            v.duration=d(rand[rank]);
            t.appendVisit(v);
            
            v.label=d(rand[rank]);
            v.duration=d(rand[rank]);
            v.duration=0;
            t.appendVisit(v);
            
            db[t.front().label].push_back(t);
            
            lSegmentLock.lock();
            localPendingSegments.push_back( std::make_pair<>(rank,db) );
            lSegmentLock.unlock();
            
            std::this_thread::sleep_for(std::chrono::microseconds(100));
            
        }
        */
    };
    
    void server(){
        
        BOOST_LOG_SEV(lg,debug) <<"#NodeManager::server";
        
        std::chrono::high_resolution_clock::time_point lastCheckpoint=std::chrono::high_resolution_clock::now();
        
        
        //send buffer for the validated batches
        std::vector<parsplice::Label> uintCommBuff2(BUFFER_SIZE);
        //receive buffer for tasks from parents and validated batches to parent
        std::vector< std::vector<parsplice::Label> > uintCommBuffv(NCOMMBUFF);
        for(int i=0;i<NCOMMBUFF;i++){
            uintCommBuffv[i]=std::vector<parsplice::Label>(BUFFER_SIZE);
        }
        

        int completed=0;
        MPI_Status status;
        
        
        //post a number of non-blocking receives
        //MPI_Irecv(&(uintCommBuff[0]),BUFFER_SIZE,MPI::UNSIGNED,MPI::ANY_SOURCE,MPI_ANY_TAG,comm,&req);
        MPI_Request reqv[NCOMMBUFF];
        for(int i=0;i<NCOMMBUFF;i++){
            MPI_Irecv(&(uintCommBuffv[i][0]),BUFFER_SIZE,MPI::UNSIGNED,MPI::ANY_SOURCE,MPI_ANY_TAG,comm,&reqv[i]);
        }
        
        while(true){
            
            //BOOST_LOG_SEV(lg,info) <<"#NodeManager::server loop";
            
            bool newTasks=false;
            bool newSegments=false;
            
            for(int i=0;i<NCOMMBUFF;i++){
                completed=0;
                
                
                //test for incoming messages
                MPI_Test(&reqv[i],&completed,&status);
                if(completed){
                    //process message
                    int peer=status.MPI_SOURCE;
                    int tag=status.MPI_TAG;
                    int count=0;
                    BOOST_LOG_SEV(lg,debug) <<"#NodeManager::server Received message "<<tag<<" from "<<peer<<" in buffer "<<i;
                    MPI_Get_count( &status, MPI::UNSIGNED, &count );
                    
                    
                    if(tag == PARSPLICE_TASKS){
                        BOOST_LOG_SEV(lg,debug) <<"#NodeManager::server PARSPLICE_TASKS "<<peer;
                        //this is a set of task requests, forward it to children
                        for(auto it=children.begin();it!=children.end();it++){
                            
                            //start a non-blocking send
                            MPI_Request r;
                            MPI_Status ss;
                            MPI_Issend(&(uintCommBuffv[i][0]),count,MPI::UNSIGNED,*it,PARSPLICE_TASKS,comm, &r);
                            
                            std::chrono::high_resolution_clock::time_point start=std::chrono::high_resolution_clock::now();
                            std::chrono::milliseconds timeout=std::chrono::milliseconds(100);
                            int completed=0;
                            int cancelled=0;
                            
                            while(true){
                                if( std::chrono::high_resolution_clock::now() - start > timeout ){
                                    MPI_Cancel(&r);
                                    MPI_Wait(&r, &ss);
                                    MPI_Test_cancelled(&ss, &cancelled);
                                    if(!cancelled){
                                        completed=1;
                                    }
                                    if(cancelled){
                                        completed=0;
                                    }
                                    BOOST_LOG_SEV(lg,error) <<"#NodeManager::server cancelling tasks forwarding "<<*it;
                                }
                                else{
                                    MPI_Test(&r, &completed, &ss);
                                }
                                
                                
                                if(completed){
                                    BOOST_LOG_SEV(lg,debug) <<"#NodeManager::server tasks forwarded "<<*it;
                                    break;
                                }
                                if(cancelled){
                                    BOOST_LOG_SEV(lg,error) <<"#NodeManager::server tasks forwarding failed "<<*it;
                                    break;
                                }
                                
                            }
                            
                            
                            //MPI_Send(&(uintCommBuffv[i][0]),count,MPI::UNSIGNED,*it,PARSPLICE_TASKS,comm);
                        }
                        
                        tasks.clear();
                        weights.clear();
                        for(int ii=0;ii<count;ii+=2){
                            tasks.push_back(uintCommBuffv[i][ii]);
                            weights.push_back(uintCommBuffv[i][ii+1]);
                        }
                        newTasks=true;
                        
                        {
                            logging::record rec = lg.open_record(keywords::severity = trace);
                            if (rec)
                            {
                                logging::record_ostream strm(rec);
                                strm<<"========== TASKS AT TIME "<<std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now()-start).count()<<"\n";
                                strm<<"========== TASKS ==========\n";
                                for(int ii=0;ii<count;ii+=2){
                                    strm<<uintCommBuffv[i][ii]<<" # "<<uintCommBuffv[i][ii+1]<<"\n";
                                }
                                
                                strm<<"========== IN PREFETCH ==========\n";
                                for(auto it=inPrefetchBuffer.begin();it!=inPrefetchBuffer.end();it++){
                                    strm<<*it<<"\n";
                                    
                                }
                                strm<<"========== IN WORK ==========\n";
                                for(auto it=inWorkBuffer.begin();it!=inWorkBuffer.end();it++){
                                    strm<<*it<<"\n";
                                }
                                
                                strm.flush();
                                lg.push_record(boost::move(rec));
                            }
                            
                        }
                    }
                    
                    if(tag == PARSPLICE_SEGMENTS){
                        BOOST_LOG_SEV(lg,debug) <<"#NodeManager::server PARSPLICE_SEGMENTS "<<peer;
                        //this is a completed batch of segments
                        int ii=0;
                        parsplice::TrajDB batch;
                        batch=deserializeDB(uintCommBuffv[i], ii);
                        logTrajDBTr(lg,batch);
                        
                        globalPendingSegments.push_back( std::make_pair<>(peer,batch) );
                        newSegments=true;
                    }
                    //post the next receive
                    MPI_Irecv(&(uintCommBuffv[i][0]),BUFFER_SIZE,MPI::UNSIGNED,MPI::ANY_SOURCE,MPI_ANY_TAG,comm,&reqv[i]);
                }
            }
            
            
            //update the db
            while(sdb->serverSingle()){};
                        
            //local validation of new segments
            //workers push segments to localPendingSegments
            lSegmentLock.lock();
            while(localPendingSegments.size()>0){
                BOOST_LOG_SEV(lg,debug) <<"#NodeManager::validate local "<<localPendingSegments.front().first;
                localValidator.validate(localPendingSegments.front().first, localPendingSegments.front().second);
                localPendingSegments.pop_front();
            }
            lSegmentLock.unlock();
            
            //release local segments to global validator
            while(true){
                parsplice::TrajDB lrelease=localValidator.release();
                if(lrelease.size()>0){
                    BOOST_LOG_SEV(lg,debug) <<"#NodeManager::validate release local ";
                    logTrajDB(lg,lrelease);
                    //release to globalValidator
                    globalValidator.validate(rank,lrelease);
                }
                else{
                    break;
                }
            }
            
            
            //push segments to the global validator
            while(globalPendingSegments.size()>0){
                BOOST_LOG_SEV(lg,debug) <<"#NodeManager::validate global "<<globalPendingSegments.front().first;
                globalValidator.validate(globalPendingSegments.front().first, globalPendingSegments.front().second);
                globalPendingSegments.pop_front();
            }
            
            //release from global validator
            while(true){
                parsplice::TrajDB release=globalValidator.release();
                if(release.size()>0){
                    BOOST_LOG_SEV(lg,debug) <<"#NodeManager::validate global release batch";
                    releaseQueue.push_back(release);
                }
                else{
                    break;
                }
            }
            
            //release segments to parent
            while(releaseQueue.size()>0){
                BOOST_LOG_SEV(lg,debug) <<"#NodeManager::server validate global release to parent "<<parent;
                
                //try sending a validated batch to the parent
                int i=0;
                serializeDB(releaseQueue.front(), uintCommBuff2, i);
                if(parent>=0){
                    
                    //start a non-blocking send
                    MPI_Request r;
                    MPI_Status ss;
                    MPI_Issend(&(uintCommBuff2[0]),i,MPI::UNSIGNED,parent,PARSPLICE_SEGMENTS,comm, &r);
                    
                    std::chrono::high_resolution_clock::time_point start=std::chrono::high_resolution_clock::now();
                    std::chrono::milliseconds timeout=std::chrono::milliseconds(100);
                    int completed=0;
                    int cancelled=0;
                    
                    while(true){
                        if( std::chrono::high_resolution_clock::now() - start > timeout ){
                            MPI_Cancel(&r);
                            MPI_Wait(&r, &ss);
                            MPI_Test_cancelled(&ss, &cancelled);
                            if(!cancelled){
                                completed=1;
                            }
                            if(cancelled){
                                completed=0;
                            }
                            BOOST_LOG_SEV(lg,error) <<"#NodeManager::server validate cancelling release to parent ";
                        }
                        else{
                            MPI_Test(&r, &completed, &ss);
                        }
                        
                        if(completed){
                            BOOST_LOG_SEV(lg,debug) <<"#NodeManager::server validate global released ";
                            break;
                        }
                        if(cancelled){
                            BOOST_LOG_SEV(lg,error) <<"#NodeManager::server validate global release failed ";
                            break;
                        }
                        
                    }
                    //this is done
                    if(completed){
                        releaseQueue.pop_front();
                    }
                    //try again later
                    if(cancelled){
                        break;
                    }
                    
                    //MPI_Send(&(uintCommBuff2[0]),i,MPI::UNSIGNED,parent,PARSPLICE_SEGMENTS,comm);
                }
                
                
            }
            
            
            /*
            while(true){
                parsplice::TrajDB release=globalValidator.release();
                if(release.size()>0){
                    BOOST_LOG_SEV(lg,debug) <<"#NodeManager::validate global release";
                    logTrajDB(lg,release);
                    //release to parent
                    int i=0;
                    serializeDB(release, uintCommBuff2, i);
                    if(parent>=0){
                        MPI_Send(&(uintCommBuff2[0]),i,MPI::UNSIGNED,parent,PARSPLICE_SEGMENTS,comm);
                    }
                    BOOST_LOG_SEV(lg,debug) <<"#NodeManager::validate global release done";
                }
                else{
                    break;
                }
            }
            */
        
            //promote fulfilled prefetches to work status
            std::set<int> promoted;
            for(int i=0;i<prefetchBuffer.size();i++){
                //is the prefetch ready?
                if(prefetchBuffer[i].second.wait_for(std::chrono::seconds(0)) == std::future_status::ready){
                    Rd data=prefetchBuffer[i].second.get();
                    workLock.lock();
                    workBuffer.push_back( std::make_pair<>(prefetchBuffer[i].first,data) );
                    inWorkBuffer.insert(prefetchBuffer[i].first);
                    BOOST_LOG_SEV(lg,info) <<"#NodeManager::server promoting state "<<prefetchBuffer[i].first;
                    if(workBuffer.size()>nWork){
                        inWorkBuffer.erase(workBuffer.front().first);
                        workBuffer.pop_front();
                    }
                    workLock.unlock();
                    promoted.insert(i);
                }
            }
            for(auto it=promoted.rbegin();it!=promoted.rend();it++){
                inPrefetchBuffer.erase( prefetchBuffer[*it].first );
                prefetchBuffer.erase(prefetchBuffer.begin()+*it);
            }
           
            
            if(newTasks){
                //if we created vacancies in the prefetch pool, fill them with new entries
                boost::random::discrete_distribution<int> dw(weights);
                int nFail=0;
                int nFailMax=1;
                //while(tasks.size()>0 and prefetchBuffer.size()<nPrefetch)
                
                for(int j=0;j<nPrefetch-prefetchBuffer.size();j++){
                    //BOOST_LOG_SEV(lg,info) <<"#NodeManager::server prefetching "<<prefetchBuffer.size()<<" of "<<nPrefetch;
                    unsigned key=tasks[dw(rand[0])];
                    //don't prefetch if already in prefetch or work buffers
                    if(inPrefetchBuffer.count(key)==0 and inWorkBuffer.count(key)==0 ){
                        BOOST_LOG_SEV(lg,info) <<"#NodeManager::server prefetching state "<<key<<" buffer state: "<<prefetchBuffer.size()<<" of "<<nPrefetch;
                        std::shared_future<Rd> f= sdb->get(DDS_GET | DDS_FORWARD,1, key);
                        prefetchBuffer.push_back(std::make_pair<>( key, f) );
                        inPrefetchBuffer.insert(key);
                    }
                    else{
                        //BOOST_LOG_SEV(lg,info) <<"#NodeManager::server prefetch rejected "<<key;
                        nFail++;
                        if(nFail>=nFailMax){
                            break;
                        }
                    }
                }
            }
            
            //checkpoint at given intervals
            if(std::chrono::high_resolution_clock::now() - lastCheckpoint> checkpointDelay  ){
                sdb->sync();
                hdb->sync();
                lastCheckpoint=std::chrono::high_resolution_clock::now();
            }
            
            
        }
        
        
        
    };
    
private:
    int nPrefetch;
    int nWork;
    int nWorkers;
    std::set<int> children;
    parsplice::SpinLock lSegmentLock;
    parsplice::SpinLock workLock;
    std::vector<int> tasks;
    std::vector<double> weights;
    std::list< std::pair<int, parsplice::TrajDB> > localPendingSegments;
    std::list< std::pair<int, parsplice::TrajDB> > globalPendingSegments;
    std::list< parsplice::TrajDB > releaseQueue;
    
    std::deque< std::pair<parsplice::Label, std::shared_future<Rd> > > prefetchBuffer;
    std::deque< std::pair<parsplice::Label, Rd> > workBuffer;
    
    std::set<parsplice::Label> inPrefetchBuffer;
    std::set<parsplice::Label> inWorkBuffer;
    
    parsplice::TrajDB ready;
    int rank;
    int parent;
    int batchSize;
    MPI_Comm comm;
    boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg;
    
    Validator localValidator;
    Validator globalValidator;
    //TimedValidator globalValidator;
    
    std::vector< std::thread > workers;
    boost::random::mt19937 rand[64];
    
    //AbstractDDS &dds;
    //STLLocalDataStore hdds;
    
    AbstractDDS *sdb;
    AbstractLocalDataStore *hdb;
    
    bool threadMultiple;
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::milliseconds checkpointDelay;
    
    ParSpliceWorkerDriver *drivers;
};





#endif /* defined(__ParSplice__NodeManager__) */
