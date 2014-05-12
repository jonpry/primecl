#include"JHlib.h"
#ifndef _WIN32
#include <pthread.h>
#include <cstdlib>
#endif

typedef struct  
{
#ifdef _WIN32
	CRITICAL_SECTION criticalSection;
#else
    pthread_mutex_t criticalSection;
#endif
	hashTable_t ht_msgQueues;
	uint32 uniqueNameCounter;
}messageQueueEnvironment_t;

messageQueueEnvironment_t messageQueueEnvironment;


/*
 * Initializes the global message queue environment
 * Call this before any other msgQueue_.. function
 */
void msgQueue_init()
{
#ifdef _WIN32
	InitializeCriticalSection(&messageQueueEnvironment.criticalSection);
#else
    pthread_mutex_init(&messageQueueEnvironment.criticalSection, NULL);
#endif
	hashTable_init(&messageQueueEnvironment.ht_msgQueues, 4);
	messageQueueEnvironment.uniqueNameCounter = 0x70000000;
}

/*
 * Creates an unconfigured and non-active message queue.
 */
#ifdef _WIN32
msgQueue_t* msgQueue_create(sint32 nameId, void (JHCALLBACK *messageProc)(msgQueue_t *msgQueue, sint32 msgId, uint32 param1, uint32 param2, void* data))
#else
msgQueue_t* msgQueue_create(sint32 nameId, void *messageProc(msgQueue_t *msgQueue, sint32 msgId, uint32 param1, uint32 param2, void* data))
#endif
{
	msgQueue_t *msgQueue = (msgQueue_t*)malloc(sizeof(msgQueue_t));
	memset(msgQueue, 0, sizeof(msgQueue_t));
	// setup
#ifdef _WIN32
	InitializeCriticalSection(&msgQueue->criticalSection);
#else
    pthread_mutex_init(&msgQueue->criticalSection, NULL);
#endif
	msgQueue->first = NULL;
	msgQueue->last = NULL;
	msgQueue->nameId = nameId;
	msgQueue->messageProc = messageProc;
	// auto-activate
#ifdef _WIN32
	EnterCriticalSection(&messageQueueEnvironment.criticalSection);
	hashTable_set(&messageQueueEnvironment.ht_msgQueues, msgQueue->nameId, msgQueue);
	LeaveCriticalSection(&messageQueueEnvironment.criticalSection);
#else
  pthread_mutex_lock(&messageQueueEnvironment.criticalSection);
  hashTable_set(&messageQueueEnvironment.ht_msgQueues, msgQueue->nameId, msgQueue);
  pthread_mutex_unlock(&messageQueueEnvironment.criticalSection);
#endif
	// return queue
	return msgQueue;
}


/*
 * Generates a unique nameId
 * Generated nameIds are always negative
 * Returns 0 on error
 */
sint32 msgQueue_generateUniqueNameId()
{
	//EnterCriticalSection(&messageQueueEnvironment.criticalSection);
	//for(sint32 i=-1; i>-100000; i--)
	//{
	//	if( hashTable_get(&messageQueueEnvironment.ht_msgQueues, i) == NULL )
	//	{
	//		LeaveCriticalSection(&messageQueueEnvironment.criticalSection);
	//		return i;
	//	}
	//}
	//LeaveCriticalSection(&messageQueueEnvironment.criticalSection);
	//return 0;
	uint32 name = 0;
#ifdef _WIN32
	EnterCriticalSection(&messageQueueEnvironment.criticalSection);
	name = messageQueueEnvironment.uniqueNameCounter;
	messageQueueEnvironment.uniqueNameCounter++;
	LeaveCriticalSection(&messageQueueEnvironment.criticalSection);
#else
  pthread_mutex_lock(&messageQueueEnvironment.criticalSection);
  name = messageQueueEnvironment.uniqueNameCounter;
  messageQueueEnvironment.uniqueNameCounter++;
  pthread_mutex_unlock(&messageQueueEnvironment.criticalSection);
#endif
	return name;
}

/*
 * Enables the message queue to receive messages
 */
void msgQueue_activate(msgQueue_t* msgQueue)
{
#ifdef _WIN32
	EnterCriticalSection(&messageQueueEnvironment.criticalSection);
	hashTable_set(&messageQueueEnvironment.ht_msgQueues, msgQueue->nameId, msgQueue);
	LeaveCriticalSection(&messageQueueEnvironment.criticalSection);
#else
  pthread_mutex_lock(&messageQueueEnvironment.criticalSection);
  hashTable_set(&messageQueueEnvironment.ht_msgQueues, msgQueue->nameId, msgQueue);
  pthread_mutex_unlock(&messageQueueEnvironment.criticalSection);
#endif
}

/*
 * Should enable the queue to only receive or block certain messages
 * maybe in future..
 */
//void msgQueue_registerFilter(messageQueue_t* msgQueue, sint32 msgId)
//{
//}

//bool msgQueue_check(msgQueue_t* msgQueue, msgInfo_t *msg)
//{
//	if( msg == NULL )
//		return false;
//	EnterCriticalSection(&msgQueue->criticalSection);
//	if( msgQueue->first )
//	{
//		msgInfoLink_t *next = msgQueue->first->next;
//		msg->msgId = msgQueue->first->msgInfo.msgId;
//		msg->paramA = msgQueue->first->msgInfo.paramA;
//		msg->paramB = msgQueue->first->msgInfo.paramB;
//		free(msgQueue->first);
//		if( next == NULL )
//		{
//			msgQueue->first = NULL;
//			msgQueue->last = NULL;
//		}
//		else
//		{
//			msgQueue->first = next;
//		}
//		LeaveCriticalSection(&msgQueue->criticalSection);
//		return true;
//	}
//	LeaveCriticalSection(&msgQueue->criticalSection);
//	return false;
//}

bool msgQueue_check(msgQueue_t* msgQueue)
{
#ifdef _WIN32
	EnterCriticalSection(&msgQueue->criticalSection);
#else
  pthread_mutex_lock(&msgQueue->criticalSection);
#endif
	if( msgQueue->first )
	{
		msgInfoLink_t *next = msgQueue->first->next;
		msgInfo_t msg;
		msg.msgId = msgQueue->first->msgInfo.msgId;
		msg.paramA = msgQueue->first->msgInfo.paramA;
		msg.paramB = msgQueue->first->msgInfo.paramB;
		msg.data = msgQueue->first->msgInfo.data;
		free(msgQueue->first);
		if( next == NULL )
		{
			msgQueue->first = NULL;
			msgQueue->last = NULL;
		}
		else
		{
			msgQueue->first = next;
		}
		msgQueue->messageProc(msgQueue, msg.msgId, msg.paramA, msg.paramB, msg.data);
		if( msg.data ) free(msg.data);
#ifdef _WIN32
		LeaveCriticalSection(&msgQueue->criticalSection);
#else
    pthread_mutex_unlock(&msgQueue->criticalSection);
#endif
		return true;
	}
#ifdef _WIN32
    LeaveCriticalSection(&msgQueue->criticalSection);
#else
    pthread_mutex_unlock(&msgQueue->criticalSection);
#endif
	return false;
}

#define MSGQUEUE_ALL	0x7FFFFFFF

bool msgQueue_postMessage(sint32 destNameId, sint32 msgId, uint32 param1, uint32 param2, void* data)
{
	if( destNameId == MSGQUEUE_ALL )
	{
		// send to all
	}
	else
	{
		// send to specific
		//EnterCriticalSection(&messageQueueEnvironment.criticalSection);
		msgQueue_t *msgQueue = (msgQueue_t*)hashTable_get(&messageQueueEnvironment.ht_msgQueues, destNameId);
		if( msgQueue == NULL )
		{
			if( data ) free(data);
			return false;
		}
		// allocate and setup message
		msgInfoLink_t *msg = (msgInfoLink_t*)malloc(sizeof(msgInfoLink_t));
		if( msg == NULL )
			return false;
		msg->msgInfo.msgId = msgId;
		msg->msgInfo.paramA = param1;
		msg->msgInfo.paramB = param2;
		msg->msgInfo.data = data;
		// append
#ifdef _WIN32
		EnterCriticalSection(&msgQueue->criticalSection);
#else
    pthread_mutex_lock(&msgQueue->criticalSection);
#endif
		if( msgQueue->last == NULL )
		{
			if( msgQueue->first != NULL )
#ifdef _WIN32
			__debugbreak(); //BUG!
#else
			raise(SIGTRAP);
#endif
			// new entry
			msg->next = NULL;
			msgQueue->first = msg;
			msgQueue->last = msg;
		}
		else
		{
			// append to last
			if( msgQueue->last->next )
#ifdef _WIN32
				__debugbreak();
#else
			    raise(SIGTRAP);
#endif
			msg->next = NULL;
			msgQueue->last->next = msg;
			msgQueue->last = msg;
		}
#ifdef _WIN32
		LeaveCriticalSection(&msgQueue->criticalSection);
#else
    pthread_mutex_unlock(&msgQueue->criticalSection);
#endif
	}
	return true;
}