#ifdef _WIN32
typedef unsigned __int32 uint32_t;
#else
#include <stdint.h>
#endif

typedef struct  
{
	int itemIndex;
}_HashTable_uint32Iterable_entry_t;

typedef struct
{
	_HashTable_uint32Iterable_entry_t *entrys;
	unsigned int *itemKeyArray;
	void **itemValueArray;
	uint32_t size;
	uint32_t count;
}hashTable_t;

void hashTable_init(hashTable_t *hashTable, int itemLimit);
void hashTable_destroy(hashTable_t *hashTable);
void hashTable_clear(hashTable_t *hashTable);
bool hashTable_set(hashTable_t *hashTable, unsigned int key, void *item);
void *hashTable_get(hashTable_t *hashTable, unsigned int key);

void** hashTable_getValueArray(hashTable_t *hashTable);
unsigned int* hashTable_getKeyArray(hashTable_t *hashTable);
unsigned int hashTable_getCount(hashTable_t *hashTable);

bool hashTable_set(hashTable_t *hashTable, char *key, void *item);
void *hashTable_get(hashTable_t *hashTable, char *key);

