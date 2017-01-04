//
//  Dispatcher.h
//  ParSplice
//
//  Created by Danny Perez on 11/11/15.
//  Copyright (c) 2015 dp. All rights reserved.
//

#ifndef __ParSplice__Dispatcher__
#define __ParSplice__Dispatcher__

#include <stdio.h>
#include <mpi.h>
#include "ParSpliceCommon.h"


class Dispatcher{
public:
    Dispatcher(MPI_Comm commp_,MPI_Comm commc_, int nTask_, int parent){
        commp=commp_;
        commc=commc_;
        nTask=nTask_;
        root=parent;
        
        
        pthread_mutex_init(&lock, NULL);
        
        pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
       
        //spin the two threads
        int ret=pthread_create(&serverTh,&attr,&Dispatcher::serverThHandle,this);
        
        ret=pthread_create(&dispatcherTh,&attr,&Dispatcher::dispatcherThHandle,this);
    };
    
    static void *serverThHandle(void *context){
        Dispatcher *p=static_cast< Dispatcher* >(context);
        p->server();
        return NULL;
    };
    
    static void *dispatcherThHandle(void *context){
        Dispatcher *p=static_cast< Dispatcher* >(context);
        p->dispatcher();
        return NULL;
    };
    
    void server(){
        std::vector<parsplice::Label> uintCommBuff(nTask);
        while(true){
            
            MPI_Bcast(&(uintCommBuff[0]), nTask, MPI::UNSIGNED , root, commp);
            BOOST_LOG_SEV(lg,info) <<"#Received tasks from "<<root;
            
            pthread_mutex_lock(&lock);
            tasks.clear();
            for(int i=0;i<nTask;i++){
                tasks.push_back(uintCommBuff[i]);
            }
            {
                logging::record rec = lg.open_record();
                if (rec)
                {
                    logging::record_ostream strm(rec);
                    strm<<"========== TASKS ==========\n";
                    
                    for(int i=0;i<nTask;i++){
                        strm<<tasks[i]<<" ";
                    }
                    
                    
                    strm.flush();
                    lg.push_record(boost::move(rec));
                }
            }
            pthread_mutex_unlock(&lock);
            
        }
    };
    
    void dispatcher(){
        BOOST_LOG_SEV(lg,info) <<"#Dispatcher dispatcher ";
        
        
        std::vector<parsplice::Label> uintCommBuff(1);
        
        MPI_Request req;
        int completed=0;
        MPI_Status status;
        
        //post an initial non-blocking receive
        MPI_Irecv(&(uintCommBuff[0]),1,MPI::UNSIGNED,MPI::ANY_SOURCE,MPI_ANY_TAG,commc,&req);
        while(true){
            MPI_Test(&req,&completed,&status);
            
            if(completed){
                int producerID=status.MPI_SOURCE;
                requestQueue.push_back(producerID);
                BOOST_LOG_SEV(lg,info) <<"#Received request from "<<producerID;
                
                //post the next receive
                completed=0;
                MPI_Irecv(&(uintCommBuff[0]),1,MPI::UNSIGNED,MPI::ANY_SOURCE,MPI_ANY_TAG,commc,&req);
            }
            
            BOOST_LOG_SEV(lg,info) <<"#Processing "<< requestQueue.size() << " requests";
            
            //process
            pthread_mutex_lock(&lock);
            if(tasks.size()>0){
                //randomly sample with replacement
                for(std::list<int>::iterator it=requestQueue.begin();it!=requestQueue.end();it++){
                    boost::random::uniform_int_distribution<int> d(0,int(tasks.size()-1));
                    uintCommBuff[0]=tasks[d(rand)];
                    MPI_Send(&(uintCommBuff[0]),1,MPI::UNSIGNED,*it,1,commc);
                    BOOST_LOG_SEV(lg,info) <<"#Processed request from "<<*it;
                }
                requestQueue.clear();
            }
            pthread_mutex_unlock(&lock);
            
        }
        
    };
protected:
    int nTask;
    int root;
    std::deque<parsplice::Label> tasks;
    boost::random::mt19937 rand;
    std::list<int> requestQueue;
    MPI_Comm commc;
    MPI_Comm commp;
    
    
    pthread_t serverTh;
    pthread_t dispatcherTh;
    pthread_mutex_t lock;
    boost::log::sources::severity_logger< boost::log::trivial::severity_level > lg;
};
#endif /* defined(__ParSplice__Dispatcher__) */
