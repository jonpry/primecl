#include"global.h"
#include "ticker.h"
#ifndef _WIN32
#include <errno.h>
#endif

extern commandlineInput_t commandlineInput;

#ifdef _WIN32
SOCKET xptClient_openConnection(char *IP, int Port)
{
	SOCKET s=socket(AF_INET,SOCK_STREAM,IPPROTO_TCP);
	SOCKADDR_IN addr;
	memset(&addr,0,sizeof(SOCKADDR_IN));
	addr.sin_family=AF_INET;
	addr.sin_port=htons(Port);
	addr.sin_addr.s_addr=inet_addr(IP);
	int result = connect(s,(SOCKADDR*)&addr,sizeof(SOCKADDR_IN));
  if( result )
	{
		return 0;
	}
#else
int xptClient_openConnection(char *IP, int Port)
{
  int s = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
  struct sockaddr_in addr;
  memset(&addr, 0, sizeof(sockaddr_in));
  addr.sin_family = AF_INET;
  addr.sin_port=htons(Port);
  addr.sin_addr.s_addr = inet_addr(IP);
  int result = connect(s, (sockaddr*)&addr, sizeof(sockaddr_in));
#endif
  if( result < 0)
{
		return 0;
	}
	return s;
}

/*
* Sends a ping request, the server will respond with the same data as fast as possible
* To measure latency we send a high precision timestamp
*/
 void xptClient_sendPing(xptClient_t* xptClient)
 {

   uint64 timestamp = getTimeHighRes();
   // build the packet
   bool sendError = false;
   xptPacketbuffer_beginWritePacket(xptClient->sendBuffer, XPT_OPC_C_PING);
   // timestamp
   xptPacketbuffer_writeU64(xptClient->sendBuffer, &sendError, timestamp);
   // finalize
   xptPacketbuffer_finalizeWritePacket(xptClient->sendBuffer);
   // send to client
   send(xptClient->clientSocket, (const char*)(xptClient->sendBuffer->buffer), xptClient->sendBuffer->parserIndex, 0);
 }
 

/*
 * Opens a new x.pushthrough connection
 * target is the server address + worker login data to use for connecting
 */
xptClient_t* xptClient_connect(jsonRequestTarget_t* target, uint32 payloadNum)
{
	// first try to connect to the given host/port
#ifdef _WIN32
	SOCKET clientSocket = xptClient_openConnection(target->ip, target->port);
#else
  int clientSocket = xptClient_openConnection(target->ip, target->port);
#endif
	if( clientSocket == 0 )
		return NULL;
#ifdef _WIN32
	// set socket as non-blocking
	unsigned int nonblocking=1;
	unsigned int cbRet;
	WSAIoctl(clientSocket, FIONBIO, &nonblocking, sizeof(nonblocking), NULL, 0, (LPDWORD)&cbRet, NULL, NULL);
#else
  int flags, err;
  flags = fcntl(clientSocket, F_GETFL, 0); 
  flags |= O_NONBLOCK;
  err = fcntl(clientSocket, F_SETFL, flags); //ignore errors for now..
#endif
	// initialize the client object
	xptClient_t* xptClient = (xptClient_t*)malloc(sizeof(xptClient_t));
	memset(xptClient, 0x00, sizeof(xptClient_t));
	xptClient->clientSocket = clientSocket;
	xptClient->sendBuffer = xptPacketbuffer_create(64*1024);
	xptClient->recvBuffer = xptPacketbuffer_create(64*1024);
	fStrCpy(xptClient->username, target->authUser, 127);
	fStrCpy(xptClient->password, target->authPass, 127);
	xptClient->payloadNum = std::max<uint32>(1, std::min<uint32>(127, payloadNum));
#ifdef _WIN32
	InitializeCriticalSection(&xptClient->cs_shareSubmit);
#else
  pthread_mutex_init(&xptClient->cs_shareSubmit, NULL);
#endif
	xptClient->list_shareSubmitQueue = simpleList_create(4);
	// send worker login
	xptClient_sendWorkerLogin(xptClient);
	// return client object
	return xptClient;
}

/*
 * Disconnects and frees the xptClient instance
 */
void xptClient_free(xptClient_t* xptClient)
{
	xptPacketbuffer_free(xptClient->sendBuffer);
	xptPacketbuffer_free(xptClient->recvBuffer);
	if( xptClient->clientSocket != 0 )
	{
#ifdef _WIN32
		closesocket(xptClient->clientSocket);
#else
    close(xptClient->clientSocket);
#endif
	}
	simpleList_free(xptClient->list_shareSubmitQueue);
	free(xptClient);
}

const int8_t base58Decode[] =
{
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
-1, 0, 1, 2, 3, 4, 5, 6, 7, 8,-1,-1,-1,-1,-1,-1,
-1, 9,10,11,12,13,14,15,16,-1,17,18,19,20,21,-1,
22,23,24,25,26,27,28,29,30,31,32,-1,-1,-1,-1,-1,
-1,33,34,35,36,37,38,39,40,41,42,43,-1,44,45,46,
47,48,49,50,51,52,53,54,55,56,57,-1,-1,-1,-1,-1,
};

/*
* Utility function to decode base58 wallet address
* dataOut should have at least 1/2 the size of base58Input
* inputLength must not exceed 200
*/
bool xptClient_decodeBase58(const char* base58Input, sint32 inputLength, uint8* dataOut, sint32* dataOutLength)
{
if( inputLength == 0 )
return false;
if( inputLength > 200 )
return false;
sint32 writeIndex = 0;
uint32 baseArray[32];
uint32 baseTrack[32];
memset(baseArray, 0x00, sizeof(baseArray));
memset(baseTrack, 0x00, sizeof(baseTrack));
uint32 baseArraySize = 1;
baseArray[0] = 0;
baseTrack[0] = 57;
// calculate exact size of output
for(sint32 i=0; i<inputLength-1; i++)
{
// multiply baseTrack with 58
for(sint32 b=baseArraySize-1; b>=0; b--)
{
uint64 multiplyWithCarry = (uint64)baseTrack[b] * 58ULL;
baseTrack[b] = (uint32)(multiplyWithCarry&0xFFFFFFFFUL);
multiplyWithCarry >>= 32;
if( multiplyWithCarry != 0 )
{
// add carry
for(sint32 carryIndex=b+1; carryIndex<baseArraySize; carryIndex++)
{
multiplyWithCarry += (uint64)baseTrack[carryIndex];
baseTrack[carryIndex] = (uint32)(multiplyWithCarry&0xFFFFFFFFUL);
multiplyWithCarry >>= 32;
if( multiplyWithCarry == 0 )
break;
}
if( multiplyWithCarry )
{
// extend
baseTrack[baseArraySize] = (uint32)multiplyWithCarry;
baseArraySize++;
}
}
}
}
// get length of output data
sint32 outputLength = 0;
uint64 last = baseTrack[baseArraySize-1];
if( last&0xFF000000 )
outputLength = baseArraySize*4;
else if( last&0xFF0000 )
outputLength = baseArraySize*4-1;
else if( last&0xFF00 )
outputLength = baseArraySize*4-2;
else
outputLength = baseArraySize*4-3;
// convert base
for(sint32 i=0; i<inputLength; i++)
{
if( base58Input[i] >= sizeof(base58Decode)/sizeof(base58Decode[0]) )
return false;
int8_t digit = base58Decode[base58Input[i]];
if( digit == -1 )
return false;
// multiply baseArray with 58
for(sint32 b=baseArraySize-1; b>=0; b--)
{
uint64 multiplyWithCarry = (uint64)baseArray[b] * 58ULL;
baseArray[b] = (uint32)(multiplyWithCarry&0xFFFFFFFFUL);
multiplyWithCarry >>= 32;
if( multiplyWithCarry != 0 )
{
// add carry
for(sint32 carryIndex=b+1; carryIndex<baseArraySize; carryIndex++)
{
multiplyWithCarry += (uint64)baseArray[carryIndex];
baseArray[carryIndex] = (uint32)(multiplyWithCarry&0xFFFFFFFFUL);
multiplyWithCarry >>= 32;
if( multiplyWithCarry == 0 )
break;
}
if( multiplyWithCarry )
{
// extend
baseArray[baseArraySize] = (uint32)multiplyWithCarry;
baseArraySize++;
}
}
}
// add base58 digit to baseArray with carry
uint64 addWithCarry = (uint64)digit;
for(sint32 b=0; addWithCarry != 0 && b<baseArraySize; b++)
{
addWithCarry += (uint64)baseArray[b];
baseArray[b] = (uint32)(addWithCarry&0xFFFFFFFFUL);
addWithCarry >>= 32;
}
if( addWithCarry )
{
// extend
baseArray[baseArraySize] = (uint32)addWithCarry;
baseArraySize++;
}
}
*dataOutLength = outputLength;
// write bytes to about
for(sint32 i=0; i<outputLength; i++)
{
dataOut[outputLength-i-1] = (uint8)(baseArray[i>>2]>>8*(i&3));
}
return true;
}

const char *version = "JP Enhanced Miner 1.0";

/*
 * Sends the worker login packet
 */
void xptClient_sendWorkerLogin(xptClient_t* xptClient)
{
	// build the packet
	bool sendError = false;
	xptPacketbuffer_beginWritePacket(xptClient->sendBuffer, XPT_OPC_C_AUTH_REQ);
	if(commandlineInput.xpt6)
		xptPacketbuffer_writeU32(xptClient->sendBuffer, &sendError, 6);								// version
	else
		xptPacketbuffer_writeU32(xptClient->sendBuffer, &sendError, 5);								// version
	
	xptPacketbuffer_writeString(xptClient->sendBuffer, xptClient->username, 128, &sendError);	// username
	xptPacketbuffer_writeString(xptClient->sendBuffer, xptClient->password, 128, &sendError);	// password
	if(!commandlineInput.xpt6)
		xptPacketbuffer_writeU32(xptClient->sendBuffer, &sendError, xptClient->payloadNum);			// payloadNum
	// write worker version to server
	xptPacketbuffer_writeString(xptClient->sendBuffer, (char*)version, 45, &sendError);	// minerVersionString
	// finalize

	uint8 walletAddressRaw[256];
	sint32 walletAddressRawLength = sizeof(walletAddressRaw);
	const char *walletAddress="AL9U76JH78egj4JA8VuiFT2x11qNYfMUFd";
	xptClient_decodeBase58(walletAddress, strlen(walletAddress), walletAddressRaw, &walletAddressRawLength);

	if(commandlineInput.xpt6){
		xptPacketbuffer_writeU8(xptClient->sendBuffer, &sendError, 1);
		xptPacketbuffer_writeU16(xptClient->sendBuffer, &sendError, 5000);
		xptPacketbuffer_writeData(xptClient->sendBuffer, walletAddressRaw+1, 20, &sendError);
	}

	xptPacketbuffer_finalizeWritePacket(xptClient->sendBuffer);
	// send to client
	send(xptClient->clientSocket, (const char*)(xptClient->sendBuffer->buffer), xptClient->sendBuffer->parserIndex, 0);
}

/*
 * Sends the share packet
 */
void xptClient_sendShare(xptClient_t* xptClient, xptShareToSubmit_t* xptShareToSubmit)
{
	//printf("Send share\n");
	// build the packet
	bool sendError = false;
	xptPacketbuffer_beginWritePacket(xptClient->sendBuffer, XPT_OPC_C_SUBMIT_SHARE);
	xptPacketbuffer_writeData(xptClient->sendBuffer, xptShareToSubmit->merkleRoot, 32, &sendError);		// merkleRoot
	xptPacketbuffer_writeData(xptClient->sendBuffer, xptShareToSubmit->prevBlockHash, 32, &sendError);	// prevBlock
	xptPacketbuffer_writeU32(xptClient->sendBuffer, &sendError, xptShareToSubmit->version);				// version
	xptPacketbuffer_writeU32(xptClient->sendBuffer, &sendError, xptShareToSubmit->nTime);				// nTime
	xptPacketbuffer_writeU32(xptClient->sendBuffer, &sendError, xptShareToSubmit->nonce);				// nNonce
	xptPacketbuffer_writeU32(xptClient->sendBuffer, &sendError, xptShareToSubmit->nBits);				// nBits
	xptPacketbuffer_writeU32(xptClient->sendBuffer, &sendError, xptShareToSubmit->sieveSize);			// sieveSize
	xptPacketbuffer_writeU32(xptClient->sendBuffer, &sendError, xptShareToSubmit->sieveCandidate);		// sieveCandidate
	// bnFixedMultiplier
	xptPacketbuffer_writeU8(xptClient->sendBuffer, &sendError, xptShareToSubmit->fixedMultiplierSize);
	xptPacketbuffer_writeData(xptClient->sendBuffer, xptShareToSubmit->fixedMultiplier, xptShareToSubmit->fixedMultiplierSize, &sendError);
	// bnChainMultiplier
	xptPacketbuffer_writeU8(xptClient->sendBuffer, &sendError, xptShareToSubmit->chainMultiplierSize);
	xptPacketbuffer_writeData(xptClient->sendBuffer, xptShareToSubmit->chainMultiplier, xptShareToSubmit->chainMultiplierSize, &sendError);
	// share id (server sends this back in shareAck, so we can identify share response)
	xptPacketbuffer_writeU32(xptClient->sendBuffer, &sendError, 0);
	// finalize
	xptPacketbuffer_finalizeWritePacket(xptClient->sendBuffer);
	// send to client
	send(xptClient->clientSocket, (const char*)(xptClient->sendBuffer->buffer), xptClient->sendBuffer->parserIndex, 0);
}

/*
 * Processes a fully received packet
 */
bool xptClient_processPacket(xptClient_t* xptClient)
{
	// printf("Received packet with opcode %d and size %d\n", xptClient->opcode, xptClient->recvSize+4);
	if( xptClient->opcode == XPT_OPC_S_AUTH_ACK )
		return xptClient_processPacket_authResponse(xptClient);
	else if( xptClient->opcode == XPT_OPC_S_WORKDATA1 )
		return xptClient_processPacket_blockData1(xptClient);
	else if( xptClient->opcode == XPT_OPC_S_SHARE_ACK )
		return xptClient_processPacket_shareAck(xptClient);
	else if( xptClient->opcode == XPT_OPC_S_PING )
		return xptClient_processPacket_ping(xptClient);


	// unknown opcodes are accepted too, for later backward compatibility
	return true;
}

/*
 * Checks for new received packets and connection events (e.g. closed connection)
 */
bool xptClient_process(xptClient_t* xptClient)
{
	
	if( xptClient == NULL )
		return false;
	// are there shares to submit?
#ifdef _WIN32
	EnterCriticalSection(&xptClient->cs_shareSubmit);
#else
    pthread_mutex_lock(&xptClient->cs_shareSubmit);
#endif
	if( xptClient->list_shareSubmitQueue->objectCount > 0 )
	{
		for(uint32 i=0; i<xptClient->list_shareSubmitQueue->objectCount; i++)
		{
			xptShareToSubmit_t* xptShareToSubmit = (xptShareToSubmit_t*)xptClient->list_shareSubmitQueue->objects[i];
			xptClient_sendShare(xptClient, xptShareToSubmit);
			free(xptShareToSubmit);
		}
		// clear list
		xptClient->list_shareSubmitQueue->objectCount = 0;
	}
#ifdef _WIN32
	LeaveCriticalSection(&xptClient->cs_shareSubmit);
#else
  pthread_mutex_unlock(&xptClient->cs_shareSubmit);
#endif

	// check if we need to send ping
   	uint32 currentTime = (uint32)time(NULL);
	if(commandlineInput.ping!=0 && xptClient->time_sentPing!=0 && currentTime > xptClient->time_sentPing + 5 && commandlineInput.usePing){
		printf("Connection stalled!\n");
		exit(0);
	}

   	if(commandlineInput.ping!=0 && xptClient->time_sendPing != 0 && currentTime >= xptClient->time_sendPing && commandlineInput.usePing)
   	{
     		xptClient_sendPing(xptClient);
		xptClient->time_sentPing = currentTime;
     		xptClient->time_sendPing = currentTime + 240; // ping every 4 minutes
   	}

	// check for packets
	uint32 packetFullSize = 4; // the packet always has at least the size of the header
	if( xptClient->recvSize > 0 )
		packetFullSize += xptClient->recvSize;
	sint32 bytesToReceive = (sint32)(packetFullSize - xptClient->recvIndex);
	// packet buffer is always large enough at this point
	sint32 r = recv(xptClient->clientSocket, (char*)(xptClient->recvBuffer->buffer+xptClient->recvIndex), bytesToReceive, 0);
	if( r < 0 )
	{
#ifdef _WIN32
		// receive error, is it a real error or just because of non blocking sockets?
		if( WSAGetLastError() != WSAEWOULDBLOCK )
		{
			xptClient->disconnected = true;
			return false;
		}
#else
    if(errno != EAGAIN)
    {
    xptClient->disconnected = true;
    return false;
    }
#endif
		return true;
	}
	if(r==0){
		xptClient->disconnected = true;
		return false;
	}
	xptClient->recvIndex += r;
	// header just received?
	if( xptClient->recvIndex == packetFullSize && packetFullSize == 4 )
	{
		// process header
		uint32 headerVal = *(uint32*)xptClient->recvBuffer->buffer;
		uint32 opcode = (headerVal&0xFF);
		uint32 packetDataSize = (headerVal>>8)&0xFFFFFF;
		// validate header size
		if( packetDataSize >= (1024*1024*2-4) )
		{
			// packets larger than 2mb are not allowed
				std::cout << "xptServer_receiveData(): Packet exceeds 2mb size limit" << std::endl;
			return false;
		}
		xptClient->recvSize = packetDataSize;
		xptClient->opcode = opcode;
		// enlarge packetBuffer if too small
		if( (xptClient->recvSize+4) > xptClient->recvBuffer->bufferLimit )
		{
			xptPacketbuffer_changeSizeLimit(xptClient->recvBuffer, (xptClient->recvSize+4));
		}
	}
	// have we received the full packet?
	if( xptClient->recvIndex >= (xptClient->recvSize+4) )
	{
		// process packet
		xptClient->recvBuffer->bufferSize = (xptClient->recvSize+4);
		if( xptClient_processPacket(xptClient) == false )
		{
			xptClient->recvIndex = 0;
			xptClient->recvSize = 0;
			xptClient->opcode = 0;
			// disconnect
			if( xptClient->clientSocket != 0 )
			{
#ifdef _WIN32
				closesocket(xptClient->clientSocket);
#else
	        close(xptClient->clientSocket);
#endif
				xptClient->clientSocket = 0;
			}
			xptClient->disconnected = true;
			return false;
		}
		xptClient->recvIndex = 0;
		xptClient->recvSize = 0;
		xptClient->opcode = 0;
	}
	// return
	return true;
}

/*
 * Returns true if the xptClient connection was disconnected from the server or should disconnect because login was invalid or awkward data received
 */
bool xptClient_isDisconnected(xptClient_t* xptClient, char** reason)
{
	if( reason )
		*reason = xptClient->disconnectReason;
	return xptClient->disconnected;
}

/*
 * Returns true if the worker login was successful
 */
bool xptClient_isAuthenticated(xptClient_t* xptClient)
{
	return (xptClient->clientState == XPT_CLIENT_STATE_LOGGED_IN);
}

void xptClient_foundShare(xptClient_t* xptClient, xptShareToSubmit_t* xptShareToSubmit)
{
#ifdef _WIN32
	EnterCriticalSection(&xptClient->cs_shareSubmit);
	simpleList_add(xptClient->list_shareSubmitQueue, xptShareToSubmit);
	LeaveCriticalSection(&xptClient->cs_shareSubmit);
#else
  pthread_mutex_lock(&xptClient->cs_shareSubmit);
  simpleList_add(xptClient->list_shareSubmitQueue, xptShareToSubmit);
  pthread_mutex_unlock(&xptClient->cs_shareSubmit);
#endif
}
