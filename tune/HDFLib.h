#include <string>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <map>
//#include <boost/filesystem.hpp>
#include <tuple>
#include <chrono>
#include <memory>

#include "hdf5.h"

herr_t file_info(hid_t loc_id, const char *name, void *opdata)
{
    H5G_stat_t statbuf;

    /*
     * Get type of the object and display its name and type.
     * The name of the object is passed to this function by 
     * the Library. Some magic :-)
     */
    H5Gget_objinfo(loc_id, name, 0, &statbuf);
    switch (statbuf.type) {
    case H5G_GROUP: 
         printf(" Object with name %s is a group \n", name);
         break;
    case H5G_DATASET: 
         printf(" Object with name %s is a dataset \n", name);
         break;
    case H5G_TYPE: 
         printf(" Object with name %s is a named datatype \n", name);
         break;
    default:
         printf(" Unable to identify an object ");
    }
    return 0;
 }

#include <fstream>

inline bool exists_test0 (const std::string& name) {
    std::ifstream f(name.c_str());
    bool good = f.good();
    f.close(); 
    return good;
}

namespace HDFLib {

#define MAX_NAME 255
//flags
const int CREATE = 1;

template<class T>
bool powerOf2(T value ) {
	if (value <= 0) {
		return false;
	}
	return (value & (value - 1)) == 0;
}


class HDFFile {
public:
	HDFFile(const std::string& filename ): _filename(filename) {
		_exists = exists_test0(filename);
        if(!_exists) {
            std::cout << "File doesn't exist: " << _exists << std::endl;    
            exit(1);
        }
        std::cout << "Exists: " << _exists << std::endl;
	}

	~HDFFile() {
		close();
	}

	void open(int flags=0) {
		//just open if file exists
		if (_exists) {

			openHDF();
			_open = true;
		}
		//else we create it
		else if (flags & CREATE) {
			createHDF();
			_open = true;
		}
		else {
			throw std::runtime_error("File does not exists and not told to create a new file");
		}
	}

	float getCompressionRatio(){
		if (!_open ||_turns==0 || _bunches==0 ) {
			return 0.0f;
		}
		return  static_cast<float>(H5Dget_storage_size(_dataset_id))/static_cast<float>(_turns*_bunches*sizeof(int16_t));
	}

	void setAttributesEnabled(bool enabled){
		if(_open){
			throw std::runtime_error("Can't change attributesEnabled after a file is open");
		}
		_attributes_enabled=enabled;

	}
	void setMinimalAttributes(bool enabled){
		if(_open){
			throw std::runtime_error("Can't change minimalAttributes after a file is open");
		}
		_minimal_attributes=enabled;

	}

	float getWriteSpeed(){
		if(_time_to_write_data==0){
			return 0.0f;
		}
		else{
			
			float megabytes=static_cast<float>( _turns*_bunches*sizeof(int16_t))/1000000.0f;
			float time_in_s=static_cast<float>(_time_to_write_data)/1000000000.0f;
			return megabytes/time_in_s;
		}
	}

	void close() {
		if (!_open) {
			return;
		}
		//close all attributes
		for (auto it : _attributes) {
			if (getAttribute(it.second) < 0) {
				continue;
			}
			_status = H5Dclose(getAttribute(it.second));
			if (_status < 0)
			{
				std::ostringstream temp;
				temp << "Attribute close failed";
				throw std::runtime_error(temp.str());
			}
			setAttribute(it.second, -1);
			_attributes [it.first] = it.second;

		}

		_status = H5Dclose(_dataset_id);
		if (_status < 0)
		{
			_status = H5Fclose(_file_id);
			std::ostringstream temp;
			temp << "Dataset close failed";
			throw std::runtime_error(temp.str());
		}

		_status = H5Pclose(_plist_id);
		if (_status < 0)
		{
			std::ostringstream temp;
			temp << "Property list close failed";
			throw std::runtime_error(temp.str());
		}

		_status = H5Sclose(_dataspace_id_data);
		if (_status < 0)
		{
			std::ostringstream temp;
			temp << "Dataspace close failed";
			throw std::runtime_error(temp.str());
		}
		if(_group_attr_id>0){
			_status = H5Gclose(_group_attr_id);
			if (_status < 0)
			{
				std::ostringstream temp;
				temp << "Group close failed";
				throw std::runtime_error(temp.str());
			}
		}

		_status = H5Gclose(_group_id);
		if (_status < 0)
		{
			std::ostringstream temp;
			temp << "Group close failed";
			throw std::runtime_error(temp.str());
		}
		if (_root_group_id >= 0) {
			_status = H5Gclose(_root_group_id);
			if (_status < 0)
			{
				std::ostringstream temp;
				temp << "Closing root group failed";
				throw std::runtime_error(temp.str());
			}
		}

		_status = H5Pclose(_file_access_plist_id);
		if(_status< 0 ){
			std::ostringstream temp;
			temp << "Closing file access property list failed";
			throw std::runtime_error(temp.str());
		}

		_status = H5Fclose(_file_id);
		if (_status < 0)
		{
			std::ostringstream temp;
			temp << "File close failed";
			throw std::runtime_error(temp.str());
		}
		H5garbage_collect();
		_open = false;
	}

	void setPlane(const std::string& plane) {
		if (_exists) {
			throw std::runtime_error("Can't change plane on a existing file");
		}
		if (_open) {
			throw std::runtime_error("Can't change plane on a open file");
		}
		if (plane != "horizontal" && plane != "vertical") {
			throw std::runtime_error("plane can only be vertical or horizontal");
		}
		_plane = plane;
		_loc = _group + "/" + _plane;
	}

	std::string getPlane(){
		return _plane;
	}

	void setBeam(const std::string& beam) {
		if (_exists) {
			throw std::runtime_error("Can't change beam on a existing file");
		}
		if (_open) {
			throw std::runtime_error("Can't change beam on a open file");
		}

		if (beam != "B1" && beam != "B2") {
			throw std::runtime_error("beam can only be B1 or B2");
		}
		_beam = beam;
		_group = "/" + _beam;
		_loc = _group + "/" + _plane;
	}

	std::string getBeam(){
		return _beam;
	}

	std::string getGroup(){
		return _group;
	}

	std::string getLoc(){
		return _loc;
	}

	void setTurns(const std::size_t& turns) {
		if (_exists) {
			throw std::runtime_error("Can't change turns on a existing file");
		}
		if (_open) {
			throw std::runtime_error("Can't change turns on a open file");
		}
		if (!powerOf2(turns)) {
			throw std::runtime_error("Turns must be a power of 2");
		}
		_setTurns(turns);

	}

	void setCompressionChunks(const std::size_t& chunkSize){
		if(chunkSize>_turns){
			_cdims[0]=_turns;
		}
		else{
			_cdims[0] = chunkSize;
		}
	}

	std::size_t getCompressionChunks(){
		return _cdims[0];
	}

	void setBunches(const std::size_t& bunches) {
		if (_exists) {
			throw std::runtime_error("Can't change bunches on a existing file");
		}
		if (_open) {
			throw std::runtime_error("Can't change bunches on a open file");
		}
		_setBunches(bunches);
	}

	std::size_t getRows() {
		if (_transpose) {
			return _bunches;
		}
		else {
			return _turns;
		}
	}

	std::size_t getColumns() {
		if (_transpose) {
			return _turns;
		}
		else {
			return _bunches;
		}
	}

	std::unique_ptr<int16_t> getData() {
		if (!_open) {
			throw std::runtime_error("File not open while trying to read data");
		}
		std::size_t size = _bunches * _turns * H5Tget_size( H5T_NATIVE_SHORT);
		void* temp_ptr = malloc(size);
		if (temp_ptr == NULL) {
			std::ostringstream temp;
			temp << "Allocating " << size << " bytes for data failed";
			throw std::runtime_error(temp.str());
		}
		_status = H5Dread(_dataset_id, H5T_NATIVE_SHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_ptr);
		if (_status < 0)
		{
			free(temp_ptr);
			std::ostringstream temp;
			temp << "Read data failed";
			throw std::runtime_error(temp.str());
		}
		return std::unique_ptr<int16_t>(reinterpret_cast<int16_t*> (temp_ptr) );
	}



	bool setData(const int16_t* data) {
		if (!_open) {
			throw std::runtime_error("File not open while trying to read data");
		}
		std::size_t start_time=getCurrentTime();
		_status = H5Dwrite(_dataset_id, H5T_NATIVE_SHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);		
		_time_to_write_data=getCurrentTime()-start_time;
		if (_status < 0)
		{
			std::ostringstream temp;
			temp << "write data failed";
			throw std::runtime_error(temp.str());
		}
		return true;
	}

	void setRowData(const int16_t* rowdata, std::size_t index) {
		_setStride(rowdata,1, _bunches, index, 0);
	}

	//will probably never be needed
	void setColumnData(const int16_t* columndata, std::size_t index) {
		_setStride(columndata,_turns, 1, 0, index);
	}

	std::unique_ptr<int16_t> operator[](int index) {
		if (_transpose) {
			return getColumnData(index);
		}
		else {
			return getRowData(index);
		}
	}

	std::unique_ptr<int16_t> getRowData(std::size_t index) {
		return _getStride(1, _bunches, index, 0);
	}


	std::unique_ptr<int16_t> getColumnData(std::size_t index) {
		return _getStride(_turns, 1, 0, index);
	}


	void setTranspose(bool value) {
		_transpose = value;
	}

	void setDefaultHeader(uint32_t header) {
		_defaultHeader = header;
	}

	void setDefaultSourceID(uint32_t sourceID) {
		_defaultSourceID = sourceID;
	}

	void setDefaultBlockSize(uint32_t BlockSize) {
		_defaultBlockSize = BlockSize;
	}

	void setDefaultTagBit(uint32_t BlockSize) {
		_defaultTagBits = BlockSize;
	}

	void setDefaultCycleNumber(uint32_t BlockSize) {
		_defaultCycleNumber = BlockSize;
	}

	void setDefaultDataWordSize(uint16_t BlockSize) {
		_defaultDataWordSize = BlockSize;
	}

	void setDefaultDataPerBunch(uint16_t BlockSize) {
		_defaultDataPerBunch = BlockSize;
	}

	void setDefaultReserved(uint32_t BlockSize) {
		_defaultReserved = BlockSize;
	}

	//headers
	std::unique_ptr<uint32_t> getHeaders() {
		return _readAttributeArray("Headers");
	}

	bool setHeaders(const uint32_t* headers) {
		auto to_return = _writeAttributeArray("Headers", headers);

		bool allOK = true;
		bool* temp = reinterpret_cast<bool*>(malloc(_turns * sizeof(bool)));
		for (std::size_t i = 0; i < _turns; i++) {
			if (headers[i] == _defaultHeader) {
				temp[i] = true;
			}
			else {
				temp[i] = false;
				allOK = false;
			}
		}
		setHeadersOK(temp);
		setHeaderOK(allOK);
		free(temp);

		return to_return;
	}

	std::unique_ptr<bool> getHeadersOK() {
		return _readAttributeArrayBool("HeadersOK");
	}

	bool setHeadersOK(const bool* headersOK) {
		return _writeAttributeArrayBool("HeadersOK", headersOK);
	}

	bool getHeaderOK() {
		return *_readAttributeArrayBool("HeaderOK");
	}

	bool setHeaderOK(const bool headerOK) {
		return _writeAttributeArrayBool("HeaderOK", &headerOK);
	}

	



//sourceID
std::unique_ptr<uint32_t> getSourceIDs() {
	return _readAttributeArray("SourceIDs");
}

bool setSourceIDs(const uint32_t* SourceIDs) {
	auto to_return = _writeAttributeArray("SourceIDs", SourceIDs);

	bool allOK = true;
	bool* temp = reinterpret_cast<bool*>(malloc(_turns * sizeof(bool)));
	for (std::size_t i = 0; i < _turns; i++) {
		if (SourceIDs[i] == _defaultSourceID) {
			temp[i] = true;
		}
		else {
			temp[i] = false;
			allOK = false;
		}
	}
	setSourceIDsOK(temp);
	setSourceIDOK(allOK);
	free(temp);
	return to_return;

}

std::unique_ptr<bool> getSourceIDsOK() {
	return _readAttributeArrayBool("SourceIDsOK");
}

bool setSourceIDsOK(const bool* SourceIDsOK) {
	return _writeAttributeArrayBool("SourceIDsOK", SourceIDsOK);
}

bool getSourceIDOK() {
	return *_readAttributeArrayBool("SourceIDOK");
}

bool setSourceIDOK(const bool SourceIDOK) {
	return _writeAttributeArrayBool("SourceIDOK", &SourceIDOK);
}

//blockSizes
std::unique_ptr<uint32_t> getBlockSizes() {
	return _readAttributeArray("BlockSizes");
}

bool setBlockSizes(const uint32_t* BlockSizes) {
	auto to_return = _writeAttributeArray("BlockSizes", BlockSizes);

	bool allOK = true;
	bool* temp = reinterpret_cast<bool*>(malloc(_turns * sizeof(bool)));
	for (std::size_t i = 0; i < _turns; i++) {
		if (BlockSizes[i] == _defaultBlockSize) {
			temp[i] = true;
		}
		else {
			temp[i] = false;
			allOK = false;
		}
	}
	setBlockSizesOK(temp);
	setBlockSizeOK(allOK);
	free(temp);
	return to_return;
}

std::unique_ptr<bool> getBlockSizesOK() {
	return _readAttributeArrayBool("BlockSizesOK");
}

bool setBlockSizesOK(const bool* BlockSizesOK) {
	return _writeAttributeArrayBool("BlockSizesOK", BlockSizesOK);
}

bool getBlockSizeOK() {
	return *_readAttributeArrayBool("BlockSizeOK");
}

bool setBlockSizeOK(const bool BlockSizeOK) {
	return _writeAttributeArrayBool("BlockSizeOK", &BlockSizeOK);
}

//TurnCounter
std::unique_ptr<uint32_t> getTurnCounters() {
	return _readAttributeArray("TurnCounters");
}

bool setTurnCounters(const uint32_t* TurnCounters) {
	auto to_return = _writeAttributeArray("TurnCounters", TurnCounters);
	bool allOK = true;
	bool* temp = reinterpret_cast<bool*>(malloc(_turns * sizeof(bool)));
	temp[_turns - 1] = true;
	for (std::size_t i = 0; i < _turns - 1; i++) {
		if (TurnCounters[i] + 1 == TurnCounters[i + 1]) {
			temp[i] = true;
		}
		else {
			temp[i] = false;
			allOK = false;
		}
	}
	setTurnCountersOK(temp);
	setTurnCounterOK(allOK);
	free(temp);
	return to_return;
}

std::unique_ptr<bool> getTurnCountersOK() {
	return _readAttributeArrayBool("TurnCountersOK");
}

bool setTurnCountersOK(const bool* TurnCountersOK) {
	return _writeAttributeArrayBool("TurnCountersOK", TurnCountersOK);
}

bool getTurnCounterOK() {
	return *_readAttributeArrayBool("TurnCounterOK");
}

bool setTurnCounterOK(const bool TurnCounterOK) {
	return _writeAttributeArrayBool("TurnCounterOK", &TurnCounterOK);
}

//TagBits
std::unique_ptr<uint32_t> getTagBits() {
	return _readAttributeArray("TagBits");
}

bool setTagBits(const uint32_t* TagBits) {
	auto to_return = _writeAttributeArray("TagBits", TagBits);
	bool allOK = true;
	bool* temp = reinterpret_cast<bool*>(malloc(_turns * sizeof(bool)));
	for (std::size_t i = 0; i < _turns; i++) {
		if (TagBits[i] == _defaultTagBits) {
			temp[i] = true;
		}
		else {
			temp[i] = false;
			allOK = false;
		}
	}
	setTagBitsOK(temp);
	setTagBitOK(allOK);
	free(temp);
	return to_return;
}

std::unique_ptr<bool> getTagBitsOK() {
	return _readAttributeArrayBool("TagBitsOK");
}

bool setTagBitsOK(const bool* TagBits) {
	return _writeAttributeArrayBool("TagBitsOK", TagBits);
}

bool getTagBitOK() {
	return *_readAttributeArrayBool("TagBitOK");
}

bool setTagBitOK(const bool TagBitOK) {
	return _writeAttributeArrayBool("TagBitOK", &TagBitOK);
}

//CycleNumbers
std::unique_ptr<uint32_t> getCycleNumbers() {
	return _readAttributeArray("CycleNumbers");
}

bool setCycleNumbers(const uint32_t* CycleNumbers) {
	auto to_return = _writeAttributeArray("CycleNumbers", CycleNumbers);

	bool allOK = true;
	bool* temp = reinterpret_cast<bool*>(malloc(_turns * sizeof(bool)));
	for (std::size_t i = 0; i < _turns; i++) {
		if (CycleNumbers[i] == _defaultCycleNumber) {
			temp[i] = true;
		}
		else {
			temp[i] = false;
			allOK = false;
		}
	}
	setCycleNumbersOK(temp);
	setCycleNumberOK(allOK);
	free(temp);
	return to_return;
}

std::unique_ptr<bool> getCycleNumbersOK() {
	return _readAttributeArrayBool("CycleNumbersOK");

}

bool setCycleNumbersOK(const bool* CycleNumbers) {
	return _writeAttributeArrayBool("CycleNumbersOK", CycleNumbers);
}

bool getCycleNumberOK() {
	return *_readAttributeArrayBool("CycleNumberOK");
}

bool setCycleNumberOK(const bool CycleNumberOK) {
	return _writeAttributeArrayBool("CycleNumberOK", &CycleNumberOK);
}


//DataWordSize
std::unique_ptr<uint16_t> getDataWordSizes() {

	return std::unique_ptr<uint16_t>( reinterpret_cast<uint16_t*> (_readAttributeArray("DataWordSizes").release()));
}

bool setDataWordSizes(const uint16_t* DataWordSizes) {
	auto to_return = _writeAttributeArray("DataWordSizes", DataWordSizes);
	bool allOK = true;
	bool* temp = reinterpret_cast<bool*>(malloc(_turns * sizeof(bool)));
	for (std::size_t i = 0; i < _turns; i++) {
		if (DataWordSizes[i] == _defaultDataWordSize) {
			temp[i] = true;
		}
		else {
			temp[i] = false;
			allOK = false;
		}
	}
	setDataWordSizesOK(temp);
	setDataWordSizeOK(allOK);
	free(temp);
	return to_return;
}

std::unique_ptr<bool> getDataWordSizesOK() {
	return _readAttributeArrayBool("DataWordSizesOK");
}

bool setDataWordSizesOK(const bool* DataWordSizes) {
	return _writeAttributeArrayBool("DataWordSizesOK", DataWordSizes);
}

bool getDataWordSizeOK() {
	return *_readAttributeArrayBool("DataWordSizeOK");
}

bool setDataWordSizeOK(const bool DataWordSizeOK) {
	return _writeAttributeArrayBool("DataWordSizeOK", &DataWordSizeOK);
}

//DataPerBunch
std::unique_ptr<uint16_t> getDataPerBunches() {

	return std::unique_ptr<uint16_t>( reinterpret_cast<uint16_t*> (_readAttributeArray("DataPerBunches").release()));
}

bool setDataPerBunches(const uint16_t* DataPerBunch) {
	auto to_return = _writeAttributeArray("DataPerBunches", DataPerBunch);
	bool allOK = true;
	bool* temp = reinterpret_cast<bool*>(malloc(_turns * sizeof(bool)));
	for (std::size_t i = 0; i < _turns; i++) {
		if (DataPerBunch[i] == _defaultDataPerBunch) {
			temp[i] = true;
		}
		else {
			temp[i] = false;
			allOK = false;
		}
	}
	setDataPerBunchesOK(temp);
	setDataPerBunchOK(allOK);
	free(temp);
	return to_return;
}

std::unique_ptr<bool> getDataPerBunchesOK() {
	return _readAttributeArrayBool("DataPerBunchesOK");
}

bool setDataPerBunchesOK(const bool* DataPerBunches) {
	return _writeAttributeArrayBool("DataPerBunchesOK", DataPerBunches);
}

bool getDataPerBunchOK() {
	return *_readAttributeArrayBool("DataPerBunchOK");
}

bool setDataPerBunchOK(const bool DataPerBunchOK) {
	return _writeAttributeArrayBool("DataPerBunchOK", &DataPerBunchOK);
}

//Reserves
std::unique_ptr<uint32_t> getReserves() {
	return _readAttributeArray("Reserves");
}

bool setReserves(const uint32_t* Reserves) {
	auto to_return = _writeAttributeArray("Reserves", Reserves);
	bool allOK = true;
	bool* temp = reinterpret_cast<bool*>(malloc(_turns * sizeof(bool)));
	for (std::size_t i = 0; i < _turns; i++) {
		if (Reserves[i] == _defaultReserved) {
			temp[i] = true;
		}
		else {
			temp[i] = false;
			allOK = false;
		}
	}
	setReservesOK(temp);
	setReserveOK(allOK);
	free(temp);
	return to_return;
}

std::unique_ptr<bool> getReservesOK() {
	return _readAttributeArrayBool("ReservesOK");
}

bool setReservesOK(const bool* Reserves) {
	return _writeAttributeArrayBool("ReservesOK", Reserves);
}

bool getReserveOK() {
	return *_readAttributeArrayBool("ReserveOK");
}

bool setReserveOK(const bool ReserveOK) {
	return _writeAttributeArrayBool("ReserveOK", &ReserveOK);
}


//CRC32
std::unique_ptr<uint32_t> getCRC32s() {
	return _readAttributeArray("CRC32s");
}

bool setCRC32s(const uint32_t* CRC32) {
	return _writeAttributeArray("CRC32s", CRC32);
}

std::unique_ptr<bool> getCRC32sOK() {
	return _readAttributeArrayBool("CRC32sOK");
}

bool setCRC32sOK(const bool* CRC32) {
	auto to_return = _writeAttributeArrayBool("CRC32sOK", CRC32);
	bool allOK = true;
	for (std::size_t i = 0; i < _turns; i++) {
		if (!CRC32[i]) {
			allOK = false;
		}
	}
	setCRC32OK(allOK);
	return to_return;
}

bool getCRC32OK() {
	return *_readAttributeArrayBool("CRC32OK");
}

bool setCRC32OK(const bool CRC32OK) {
	return _writeAttributeArrayBool("CRC32OK", &CRC32OK);
}

//APECTags
std::unique_ptr<bool> getAPECTagsOK() {
	return _readAttributeArrayBool("APECTagsOK");
}

bool setAPECTagsOK(const bool* APECTag) {
	auto to_return = _writeAttributeArrayBool("APECTagsOK", APECTag);
	bool allOK = true;
	for (std::size_t i = 0; i < _turns; i++) {
		if (!APECTag[i]) {
			allOK = false;
		}
	}
	setAPECTagOK(allOK);
	return to_return;
}

bool getAPECTagOK() {
	return *_readAttributeArrayBool("APECTagOK");
}

bool setAPECTagOK(const bool APECTagOK) {
	return _writeAttributeArrayBool("APECTagOK", &APECTagOK);
}

//NotInTables
std::unique_ptr<bool> getNotInTablesOK() {
	return _readAttributeArrayBool("NotInTablesOK");
}

bool setNotInTablesOK(const bool* NotInTable) {
	auto to_return = _writeAttributeArrayBool("NotInTablesOK", NotInTable);
	bool allOK = true;
	for (std::size_t i = 0; i < _turns; i++) {
		if (!NotInTable[i]) {
			allOK = false;
		}
	}
	setNotInTableOK(allOK);
	return to_return;
}

bool getNotInTableOK() {
	return *_readAttributeArrayBool("NotInTableOK");
}

bool setNotInTableOK(const bool NotInTableOK) {
	return _writeAttributeArrayBool("NotInTableOK", &NotInTableOK);
}

//Disparities
std::unique_ptr<bool> getDisparitiesOK() {
	return _readAttributeArrayBool("DisparitiesOK");
}

bool setDisparitiesOK(const bool* Disparities) {
	auto to_return = _writeAttributeArrayBool("DisparitiesOK", Disparities);
	bool allOK = true;
	for (std::size_t i = 0; i < _turns; i++) {
		if (!Disparities[i]) {
			allOK = false;
		}
	}
	setDisparityOK(allOK);
	return to_return;
}

bool getDisparityOK() {
	return *_readAttributeArrayBool("DisparityOK");
}

bool setDisparityOK(const bool Disparity) {
	return _writeAttributeArrayBool("DisparityOK", &Disparity);
}

//TriggerTime
uint64_t getTriggerTime() {

	return *reinterpret_cast<uint64_t*> (_readAttributeArray("TriggerTime").get());
}

bool setTriggerTime(const uint64_t triggertime) {
	return _writeAttributeArray("TriggerTime", &triggertime);
}

//FreezeTime
uint64_t getFreezeTime() {

	return *reinterpret_cast<uint64_t*> (_readAttributeArray("FreezeTime").get());
}

bool setFreezeTime(const uint64_t Freeze) {
	return _writeAttributeArray("FreezeTime", &Freeze);
}

//SaveTime
uint64_t getSaveTime() {

	return *reinterpret_cast<uint64_t*> (_readAttributeArray("SaveTime").get());
}

bool setSaveTime(const uint64_t Save) {
	return _writeAttributeArray("SaveTime", &Save);
}

//SaveTime
std::string getTriggerType() {

	return _readAttributeString("TriggerType");
}

bool setTriggerType(const std::string& TriggerType) {
	return _writeAttributeArray("TriggerType", TriggerType.c_str());
}

//CycleId
uint32_t getCycleId() {

	return *reinterpret_cast<uint32_t*> (_readAttributeArray("CycleId").get());
}

bool setCycleId(const uint32_t CycleId) {
	return _writeAttributeArray("CycleId", &CycleId);
}

//CycleName
std::string getCycleName() {

	return _readAttributeString("CycleName");
}

bool setCycleName(const std::string& CycleName) {
	return _writeAttributeArray("CycleName", CycleName.c_str());
}

//CycleStamp
uint64_t getCycleStamp() {

	return *reinterpret_cast<uint64_t*> (_readAttributeArray("CycleStamp").get());
}

bool setCycleStamp(const uint64_t CycleStamp) {
	return _writeAttributeArray("CycleStamp", &CycleStamp);
}

//extraCondition
std::string getExtraCondition() {

	return _readAttributeString("ExtraCondition");
}

bool setExtraCondition(const std::string& ExtraCondition) {
	return _writeAttributeArray("ExtraCondition", ExtraCondition.c_str());
}

//TimeStamp
uint64_t getTimeStamp() {

	return *reinterpret_cast<uint64_t*> (_readAttributeArray("TimeStamp").get());
}

bool setTimeStamp(const uint64_t TimeStamp) {
	return _writeAttributeArray("TimeStamp", &TimeStamp);
}

//TimingDomainName
std::string getTimingDomainName() {

	return _readAttributeString("TimingDomainName");
}

bool setTimingDomainName(const std::string& TimingDomainName) {
	return _writeAttributeArray("TimingDomainName", TimingDomainName.c_str());
}


protected:
private:

	bool _transpose = false;

	uint32_t _defaultHeader = 1331852081;
	uint32_t _defaultSourceID = 0;
	uint32_t _defaultBlockSize = 7128;
	uint32_t _defaultTagBits = 4096;
	uint32_t _defaultCycleNumber = 0;
	uint16_t _defaultDataWordSize = 0;
	uint16_t _defaultDataPerBunch = 0;
	uint32_t _defaultReserved = 0;

	std::size_t getCurrentTime(){
		std::chrono::time_point<std::chrono::high_resolution_clock> now = std::chrono::high_resolution_clock::now();
		auto duration = now.time_since_epoch();
		auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(duration);
		return nanoseconds.count();
	}


	//			attribute_id, dim, type
	typedef std::tuple <hid_t, hsize_t, hid_t> Attribute;
	hid_t getAttribute(Attribute& attr) {
		return std::get<0>(attr);
	}
	void setAttribute(Attribute& attr, hid_t id) {
		std::get<0>(attr) = id;
	}

	hsize_t getDim(Attribute& attr) {
		return std::get<1>(attr);
	}
	hid_t getType(Attribute& attr) {
		return std::get<2>(attr);
	}


std::unique_ptr<int16_t> _getStride(hsize_t count0, hsize_t count1, hsize_t offset0, hsize_t offset1,hsize_t stride0=1,hsize_t stride1=1,hsize_t block0=1,hsize_t block1=1) {
	if (!_open) {
        std::cout << "Getting stride" << std::endl;
		throw std::runtime_error("File not open while trying to read row");
	}
	hsize_t     count[2] = {count0, count1};
	hsize_t     offset[2] = {offset0, offset1};
	hsize_t     stride[2] = {stride0, stride1};
	hsize_t     block[2] = {block0, block1};

	int16_t* temp_ptr = reinterpret_cast<int16_t*>(malloc(count0 * count1 * sizeof(int16_t)));
	if (temp_ptr == NULL) {
		throw std::runtime_error("Allocating memory in getRowData Failed");
	}

	hid_t memspace_id = H5Screate_simple (2, count, NULL);
	if (memspace_id < 0) {
		throw std::runtime_error("Creating memspace in getRowData failed");
	}

	_status = H5Sselect_hyperslab (_dataspace_id_data, H5S_SELECT_SET, offset, stride, count, block);
	if (_status < 0) {
		throw std::runtime_error("Selecting hyperslab in getRowData failed");

	}
	_status = H5Dread (_dataset_id, H5T_NATIVE_SHORT, memspace_id, _dataspace_id_data, H5P_DEFAULT, temp_ptr);
	if (_status < 0) {
		throw std::runtime_error("Reading data in getRowData failed");
	}

	H5Sclose(memspace_id);

	return std::unique_ptr<int16_t>(temp_ptr);
}

void _setStride(const int16_t* data,hsize_t count0, hsize_t count1, hsize_t offset0, hsize_t offset1,hsize_t stride0=1,hsize_t stride1=1,hsize_t block0=1,hsize_t block1=1) {
	if (!_open) {
		throw std::runtime_error("File not open while trying to read row");
	}
	hsize_t     count[2] = {count0, count1};
	hsize_t     offset[2] = {offset0, offset1};
	hsize_t     stride[2] = {stride0, stride1};
	hsize_t     block[2] = {block0, block1};

	hid_t memspace_id = H5Screate_simple (2, count, NULL);
	if (memspace_id < 0) {
		throw std::runtime_error("Creating memspace in getRowData failed");
	}

	_status = H5Sselect_hyperslab (_dataspace_id_data, H5S_SELECT_SET, offset, stride, count, block);
	if (_status < 0) {
		throw std::runtime_error("Selecting hyperslab in getRowData failed");

	}
	_status = H5Dwrite (_dataset_id, H5T_NATIVE_SHORT, memspace_id, _dataspace_id_data, H5P_DEFAULT, data);
	if (_status < 0) {
		throw std::runtime_error("Reading data in getRowData failed");
	}

	H5Sclose(memspace_id);

}

std::unique_ptr<uint32_t> _readAttributeArray(const std::string& attr) {
	if(!_attributes_enabled){
		throw std::runtime_error("Can't read attribute if it is disabled");
	}
	if (!_open) {
		throw std::runtime_error("File not open while trying to read attribute");
	}
	if (_attributes.count(attr) != 1) {
		throw std::runtime_error("Attribute not found in internal map");
	}
	if (getAttribute(_attributes[attr]) <= 0) {
		std::ostringstream temp;
		temp << "attribute " << attr << " not found in file";
		throw std::runtime_error(temp.str());
	}

	std::size_t size = getDim(_attributes[attr]) * H5Tget_size( getType(_attributes[attr]));
	void* temp_ptr = malloc(size);
	if (temp_ptr == NULL) {
		std::ostringstream temp;
		temp << "Allocating " << size << " bytes for attribute " << attr << " failed";
		throw std::runtime_error(temp.str());
	}

	_status = H5Dread(getAttribute(_attributes[attr]), getType(_attributes[attr]), H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_ptr);
	if (_status < 0) {
		free(temp_ptr);
		std::ostringstream temp;
		temp << "Reading attribute " << attr << " failed";
		throw std::runtime_error(temp.str());
	}

	std::unique_ptr<uint32_t> temp_ptr_unique(reinterpret_cast<uint32_t*>(temp_ptr ) );

	return temp_ptr_unique;
}

bool _writeAttributeArray(const std::string& attr, const void* data) {
	if(!_attributes_enabled){
		throw std::runtime_error("Can't write attribute if it is disabled");
	}
	if (!_open) {
		throw std::runtime_error("File not open while trying to read attribute");
	}
	if (_attributes.count(attr) != 1) {
		throw std::runtime_error("Attribute not found in internal map");
	}
	if (getAttribute(_attributes[attr]) <= 0) {
		std::ostringstream temp;
		temp << "attribute " << attr << " not found in file";
		throw std::runtime_error(temp.str());
	}

	//std::size_t size = getDim(_attributes[attr]) * H5Tget_size( getType(_attributes[attr]));
	_status = H5Dwrite(getAttribute(_attributes[attr]), getType(_attributes[attr]),H5S_ALL, H5S_ALL, H5P_DEFAULT,  data);
	if (_status < 0) {
		std::ostringstream temp;
		temp << "Writing attribute " << attr << " failed";
		throw std::runtime_error(temp.str());
	}

	return true;
}

std::string _readAttributeString(const std::string& attr) {
	if(!_attributes_enabled){
		throw std::runtime_error("Can't read attribute if it is disabled");
	}
	if (!_open) {
		throw std::runtime_error("File not open while trying to read attribute");
	}
	if (_attributes.count(attr) != 1) {
		throw std::runtime_error("Attribute not found in internal map");
	}
	if (getAttribute(_attributes[attr]) <= 0) {
		std::ostringstream temp;
		temp << "attribute " << attr << " not found in file";
		throw std::runtime_error(temp.str());
	}

	std::size_t size = getDim(_attributes[attr]) * H5Tget_size( getType(_attributes[attr]));
	void* temp_ptr = malloc(size);
	if (temp_ptr == NULL) {
		std::ostringstream temp;
		temp << "Allocating " << size << " bytes for attribute " << attr << " failed";
		throw std::runtime_error(temp.str());
	}

	_status = H5Dread(getAttribute(_attributes[attr]), getType(_attributes[attr]), H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_ptr);
	std::string test(reinterpret_cast<const char*> (temp_ptr) );
	free(temp_ptr);
	if (_status < 0) {
		std::ostringstream temp;
		temp << "Reading attribute " << attr << " failed";
		throw std::runtime_error(temp.str());
	}


	return test;
}

bool _writeAttributeArrayBool(const std::string& attr, const bool* data) {
	if(!_attributes_enabled){
		throw std::runtime_error("Can't write attribute if it is disabled");
	}
	if (!H5Tequal(getType(_attributes[attr]), H5T_NATIVE_CHAR)) {
		std::ostringstream temp;
		temp << "Calling _writeAttributeArrayBool for attribute " << attr << " which is not of type H5T_NATIVE_CHAR";
		throw std::runtime_error(temp.str());
	}
	//should be char
	std::size_t size = getDim(_attributes[attr]) * H5Tget_size( getType(_attributes[attr]));
	uint8_t* temp_ptr = reinterpret_cast<uint8_t*>(malloc(size));
	for (std::size_t i = 0; i < getDim(_attributes[attr]); i++) {
		if (data[i]) {
			temp_ptr[i] = 1;
		}
		else {
			temp_ptr[i] = 0;
		}
	}
	bool temp_return = _writeAttributeArray(attr, data);
	free(temp_ptr);

	return temp_return;
}


std::unique_ptr<bool> _readAttributeArrayBool(const std::string& attr) {
	if(!_attributes_enabled){
		throw std::runtime_error("Can't read attribute if it is disabled");
	}
	if (!H5Tequal(getType(_attributes[attr]), H5T_NATIVE_CHAR)) {
		std::ostringstream temp;
		temp << "Calling _readAttributeArrayBool for attribute " << attr << " which is not of type H5T_NATIVE_CHAR";
		throw std::runtime_error(temp.str());
	}
	//get data as uint32
	std::unique_ptr<uint32_t> temp_ptr = _readAttributeArray(attr);
	//reinterpret to bool but it is actually a char array
	std::unique_ptr<bool> temp_ptr1(reinterpret_cast<bool*> (temp_ptr.release()) );
	//reinterpret is as a char array without giving up ownership
	uint8_t* temp_ptr2 = reinterpret_cast<uint8_t*>(temp_ptr1.get());
	for (std::size_t i = 0; i < getDim(_attributes[attr]); i++) {
		//these are pointing to the same memory but we read before we write and bool is smaller than uint8_t so we should be fine
		temp_ptr1.get()[i] = temp_ptr2[i] == 1;
	}

	return temp_ptr1;
}



void _setPlane(const std::string& plane) {
	if (plane != "horizontal" && plane != "vertical") {
		throw std::runtime_error("plane can only be vertical or horizontal");
	}
	_plane = plane;
	_loc = _group + "/" + _plane;
}

void _setBeam(const std::string& beam) {
	if (beam != "B1" && beam != "B2") {
		throw std::runtime_error("beam can only be B1 or B2");
	}
	_beam = beam;
	_group = "/" + _beam;
	_loc = _group + "/" + _plane;
}

void _setTurns(const std::size_t& turns) {

	if (!powerOf2(turns)) {
		throw std::runtime_error("Turns must be a power of 2");
	}
	_turns = turns;
	_dims[0] = _turns;
	if (_turns > 8192) {
		_cdims[0] = 8192;
	}
	else {
		_cdims[0] = _turns;
	}
}


void _setBunches(const std::size_t& bunches) {
	_bunches = bunches;
	_dims[1] = _bunches;
}

void _initialiseAttributeMap() {
	if(!_attributes_enabled){
		return;
	}

	_attributes["TurnCounters"] = std::make_tuple(-1, _turns, H5T_STD_U32LE);
	_attributes["TurnCountersOK"] = std::make_tuple(-1, _turns, H5T_NATIVE_CHAR);
	_attributes["TurnCounterOK"] = std::make_tuple(-1, 1, H5T_NATIVE_CHAR);
	if(_minimal_attributes){
		return;
	}

	_attributes["Headers"] = std::make_tuple(-1, _turns, H5T_STD_U32LE);
	_attributes["HeadersOK"] = std::make_tuple(-1, _turns, H5T_NATIVE_CHAR);
	_attributes["HeaderOK"] = std::make_tuple(-1, 1, H5T_NATIVE_CHAR);

	_attributes["SourceIDs"] = std::make_tuple(-1, _turns, H5T_STD_U32LE);
	_attributes["SourceIDsOK"] = std::make_tuple(-1, _turns, H5T_NATIVE_CHAR);
	_attributes["SourceIDOK"] = std::make_tuple(-1, 1, H5T_NATIVE_CHAR);

	_attributes["BlockSizes"] = std::make_tuple(-1, _turns, H5T_STD_U32LE);
	_attributes["BlockSizesOK"] = std::make_tuple(-1, _turns, H5T_NATIVE_CHAR);
	_attributes["BlockSizeOK"] = std::make_tuple(-1, 1, H5T_NATIVE_CHAR);

	_attributes["TagBits"] = std::make_tuple(-1, _turns, H5T_STD_U32LE);
	_attributes["TagBitsOK"] = std::make_tuple(-1, _turns, H5T_NATIVE_CHAR);
	_attributes["TagBitOK"] = std::make_tuple(-1, 1, H5T_NATIVE_CHAR);

	_attributes["CycleNumbers"] = std::make_tuple(-1, _turns, H5T_STD_U32LE);
	_attributes["CycleNumbersOK"] = std::make_tuple(-1, _turns, H5T_NATIVE_CHAR);
	_attributes["CycleNumberOK"] = std::make_tuple(-1, 1, H5T_NATIVE_CHAR);

	_attributes["DataWordSizes"] = std::make_tuple(-1, _turns, H5T_NATIVE_SHORT);
	_attributes["DataWordSizesOK"] = std::make_tuple(-1, _turns, H5T_NATIVE_CHAR);
	_attributes["DataWordSizeOK"] = std::make_tuple(-1, 1, H5T_NATIVE_CHAR);

	_attributes["DataPerBunches"] = std::make_tuple(-1, _turns, H5T_NATIVE_SHORT);
	_attributes["DataPerBunchesOK"] = std::make_tuple(-1, _turns, H5T_NATIVE_CHAR);
	_attributes["DataPerBunchOK"] = std::make_tuple(-1, 1, H5T_NATIVE_CHAR);

	_attributes["Reserves"] = std::make_tuple(-1, _turns, H5T_STD_U32LE);
	_attributes["ReservesOK"] = std::make_tuple(-1, _turns, H5T_NATIVE_CHAR);
	_attributes["ReserveOK"] = std::make_tuple(-1, 1, H5T_NATIVE_CHAR);

	_attributes["CRC32s"] = std::make_tuple(-1, _turns, H5T_STD_U32LE);
	_attributes["CRC32sOK"] = std::make_tuple(-1, _turns, H5T_NATIVE_CHAR);
	_attributes["CRC32OK"] = std::make_tuple(-1, 1, H5T_NATIVE_CHAR);

	_attributes["APECTagsOK"] = std::make_tuple(-1, _turns, H5T_NATIVE_CHAR);
	_attributes["APECTagOK"] = std::make_tuple(-1, 1, H5T_NATIVE_CHAR);

	_attributes["NotInTablesOK"] = std::make_tuple(-1, _turns, H5T_NATIVE_CHAR);
	_attributes["NotInTableOK"] = std::make_tuple(-1, 1, H5T_NATIVE_CHAR);

	_attributes["DisparitiesOK"] = std::make_tuple(-1, _turns, H5T_NATIVE_CHAR);
	_attributes["DisparityOK"] = std::make_tuple(-1, 1, H5T_NATIVE_CHAR);

	_attributes["TriggerTime"] = std::make_tuple(-1, 1, H5T_STD_U64LE);
	_attributes["FreezeTime"] = std::make_tuple(-1, 1, H5T_STD_U64LE);
	_attributes["SaveTime"] = std::make_tuple(-1, 1, H5T_STD_U64LE);

	_attributes["TriggerType"] = std::make_tuple(-1, MAX_NAME, H5T_NATIVE_CHAR);

	_attributes["CycleId"] = std::make_tuple(-1, 1, H5T_STD_U32LE);

	_attributes["CycleName"] = std::make_tuple(-1, MAX_NAME, H5T_NATIVE_CHAR);

	_attributes["CycleStamp"] = std::make_tuple(-1, 1, H5T_STD_U64LE);

	_attributes["ExtraCondition"] = std::make_tuple(-1, MAX_NAME, H5T_NATIVE_CHAR);

	_attributes["TimeStamp"] = std::make_tuple(-1, 1, H5T_STD_U64LE);

	_attributes["TimingDomainName"] = std::make_tuple(-1, MAX_NAME, H5T_NATIVE_CHAR);

}

void openHDF() {
	int otype = 0;
	ssize_t len;
	char group_name[MAX_NAME];
	char dataset_name[MAX_NAME];
	char attribute_name[MAX_NAME];
	hsize_t nobj;
	hid_t tid;
	int ndims;
	hsize_t dims[2];
	int na;
	hid_t atype;
	hid_t aspace;
	hsize_t adim;

	//open file
	_file_id = H5Fopen (_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    
    if (_file_id < 0)
	{
		std::ostringstream temp;
		temp << "opening of file " << _filename << " failed";
		throw std::runtime_error(temp.str());
	}

	//get file access property
	_file_access_plist_id = H5Fget_access_plist(_file_id);
	if (_file_access_plist_id < 0)
	{
		std::ostringstream temp;
		temp << "Could not get file access property list";
		throw std::runtime_error(temp.str());
	}

	//open root group
	_root_group_id = H5Gopen1(_file_id, "/");
    std::cout << "Root grp: " << _root_group_id << std::endl;
	if (_root_group_id < 0)
	{
        std::cout << "Failed???\n" << std::endl;
		std::ostringstream temp;
		temp << "opening of group / failed";
		throw std::runtime_error(temp.str());
	}


	//get number of items in root group
	_status = H5Gget_num_objs(_root_group_id, &nobj);
	if (_status < 0)
	{
		std::ostringstream temp;
		temp << "Getting number of objects for group / failed";
		throw std::runtime_error(temp.str());
	}
	if (nobj > 1) {
		throw std::runtime_error("The root group / should only contain one object");
	}


	//get the name, should be either B1 or B2
	len = H5Gget_objname_by_idx(_root_group_id, 0, group_name, (size_t)MAX_NAME );
	if (std::string(group_name) != "B1" && std::string(group_name) != "B2") {
		throw std::runtime_error("group /B1 or /B2 not found");
	}

	otype =  H5Gget_objtype_by_idx(_root_group_id, 0 );
	if (otype != H5G_GROUP) {
		std::ostringstream temp;
		temp << "The only member in the root group is not a group";
		throw std::runtime_error(temp.str());
	}

    printf(" Objects in the root group are:\n");
    printf("\n");

    H5Giterate(_file_id, "B2", NULL, file_info, NULL);

	_setBeam(group_name);

    std::cout << "Group name: " << group_name  << std::endl;

	//open the group
	_group_id = H5Gopen1(_root_group_id, group_name);
	if (_group_id < 0)
	{
		std::ostringstream temp;
		temp << "opening of group /" << group_name << " failed";
		throw std::runtime_error(temp.str());
	}


	//get number of items in main data group
	_status = H5Gget_num_objs(_group_id, &nobj);
	if (_status < 0)
	{
		std::ostringstream temp;
		temp << "Getting number of objects for group /" << group_name << " failed";
		throw std::runtime_error(temp.str());
	}
	if (nobj > 2) {
		std::ostringstream temp;
		temp << "The group /" << group_name << " should maximum contain 2 objects";
		throw std::runtime_error(temp.str());
		throw std::runtime_error("");
	}
	if(nobj==2){
		_attributes_enabled=true;
	}
	for(unsigned i=0;i<nobj;i++){
		//get the name, should be either horizontal or vertical
		len = H5Gget_objname_by_idx(_group_id, i, dataset_name, (size_t)MAX_NAME );
		if (std::string(dataset_name) != "vertical" && std::string(dataset_name) != "horizontal") {
			continue;
		}

		//check so it is a dataset
		otype =  H5Gget_objtype_by_idx(_group_id, i );
		if (otype != H5G_DATASET) {
			std::ostringstream temp;
			temp << "A member in the group " << group_name << " with the name horzontal or vertical must be a dataset";
			throw std::runtime_error(temp.str());
		}

	}
	if (std::string(dataset_name) != "vertical" && std::string(dataset_name) != "horizontal") {
		throw std::runtime_error("could not find vertical or horizontal dataset");
	}


	_setPlane(dataset_name);

	_dataset_id = H5Dopen1(_group_id, dataset_name);
	if (_dataset_id < 0)
	{
		std::ostringstream temp;
		temp << "Could not open dataset /" << group_name << "/" << dataset_name;
		throw std::runtime_error(temp.str());
	}
	_dataspace_id_data = H5Dget_space(_dataset_id);
	if (_dataspace_id_data < 0)
	{
		std::ostringstream temp;
		temp << "Could not get dataspace for /" << group_name << "/" << dataset_name;
		throw std::runtime_error(temp.str());
	}
	tid = H5Dget_type(_dataset_id);
	if (!H5Tequal(tid, H5T_NATIVE_SHORT)) {
		std::ostringstream temp;
		temp << "The datatype for dataset /" << group_name << "/" << dataset_name << " is not H5T_NATIVE_SHORT";
		throw std::runtime_error(temp.str());
	}
	H5Tclose(tid);

	_plist_id = H5Dget_create_plist(_dataset_id);
	if (_plist_id < 0)
	{
		std::ostringstream temp;
		temp << "Could not get property list for /" << group_name << "/" << dataset_name;
		throw std::runtime_error(temp.str());
	}

	ndims = H5Sget_simple_extent_ndims(_dataspace_id_data);
	if (ndims != 2)
	{
		std::ostringstream temp;
		temp << "Dimensions are not 2 in dataset /" << group_name << "/" << dataset_name;
		throw std::runtime_error(temp.str());
	}
	_status = H5Sget_simple_extent_dims(_dataspace_id_data, dims, NULL);
	if (_status < 0)
	{
		std::ostringstream temp;
		temp << "Failed to get dimensions from dataset /" << group_name << "/" << dataset_name;
		throw std::runtime_error(temp.str());
	}
	_setTurns(dims[0]);
	_setBunches(dims[1]);
	_initialiseAttributeMap();

	if(_attributes_enabled){
		//open the group
		_group_attr_id = H5Gopen1(_group_id, "attributes");

        std::cout << "ATTR ID: " << _group_attr_id << std::endl;
		if (_group_attr_id < 0)
		{
			std::ostringstream temp;
			temp << "opening of group /" << group_name<<"/attributes"<< " failed";
			throw std::runtime_error(temp.str());
		}


		//get number of items in attributes group
		_status = H5Gget_num_objs(_group_attr_id, &nobj);
		if (_status < 0)
		{
			std::ostringstream temp;
			temp << "Getting number of objects for group /" << group_name<<"/attributes" << " failed";
			throw std::runtime_error(temp.str());
		}

		for(unsigned i=0;i<nobj;i++){
			len = H5Gget_objname_by_idx(_group_attr_id, i, dataset_name, (size_t)MAX_NAME );
			//check so it is a dataset
			otype =  H5Gget_objtype_by_idx(_group_attr_id, i );
			if (otype != H5G_DATASET) {
				std::ostringstream temp;
				temp << "Member "<<dataset_name<<" in the attribute group is not a dataset";
				throw std::runtime_error(temp.str());
			}
			hid_t aid = H5Dopen1(_group_attr_id, dataset_name);
			if (aid < 0)
			{
				std::ostringstream temp;
				temp << "Could not open dataset /" << group_name << "/attributes/" << dataset_name;
				throw std::runtime_error(temp.str());
			}
			aspace = H5Dget_space(aid);
			if (aspace < 0)
			{
				std::ostringstream temp;
				temp << "Could not get dataspace for /" << group_name << "/" << dataset_name;
				throw std::runtime_error(temp.str());
			}
			atype = H5Dget_type(aid);
			ndims = H5Sget_simple_extent_ndims(aspace);
			if (ndims != 1)
			{
				std::ostringstream temp;
				temp << "Dimensions are not 1 in dataset /" << group_name << "/attributes/" << dataset_name;
				throw std::runtime_error(temp.str());
			}
			_status = H5Sget_simple_extent_dims(aspace, &adim, NULL);
			if (_status < 0)
			{
				std::ostringstream temp;
				temp << "Failed to get dimensions from dataset /" << group_name << "/attributes/" << dataset_name;
				throw std::runtime_error(temp.str());
			}
			H5Sclose(aspace);

			std::string attribute_name_str = dataset_name;
			//headers
			if (_attributes.count(attribute_name_str)) {
				//check dim
				if (adim != getDim(_attributes[attribute_name_str])) {
					std::ostringstream temp;
					temp << "Dimensions of attribute" << attribute_name_str << " for dataset /" << group_name << "/" << dataset_name << " does not match";
					throw std::runtime_error(temp.str());
				}
				//check type
				if (!H5Tequal(atype, getType(_attributes[attribute_name_str]))) {
					std::ostringstream temp;
					temp << "Type for attribute" << attribute_name_str << " for dataset /" << group_name << "/" << dataset_name << " is not H5T_STD_U32LE";
					throw std::runtime_error(temp.str());
				}
				H5Tclose(atype);
				setAttribute(_attributes[attribute_name_str], aid);
			}

			else {
				std::ostringstream temp;
				temp << "Attribute " << attribute_name_str << " is not defined in library";
				throw std::runtime_error(temp.str());
			}


		}




	}
	else{
		return;
	}

}


void createHDF() {

	//---------------------create file access property identifier---------------
	_file_access_plist_id = H5Pcreate(H5P_FILE_ACCESS);
	if (_file_access_plist_id < 0)
	{
		std::ostringstream temp;
		temp << "H5Pcreate(H5P_FILE_ACCESS) failed";
		throw std::runtime_error(temp.str());
	}
	#if  H5_VERS_MINOR == 12
		//---------------------use the latest file format---------------------------
		_status=H5Pset_libver_bounds(_file_access_plist_id, H5F_LIBVER_V18, H5F_LIBVER_V18);
	#else
		_status=H5Pset_libver_bounds(_file_access_plist_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
	#endif
	if (_status < 0)
	{
		std::ostringstream temp;
		temp << "H5Pset_attr_phase_change failed";
		throw std::runtime_error(temp.str());
	}

	//---------------------create file------------------------------------------
	_file_id = H5Fcreate(_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, _file_access_plist_id);
	if (_file_id < 0)
	{
		std::ostringstream temp;
		temp << "creation of file " << _filename << " failed";
		throw std::runtime_error(temp.str());
	}

	//-----------------------create group----------------------------------------

	_group_id = H5Gcreate2(_file_id, _group.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (_group_id < 0)
	{
		std::ostringstream temp;
		temp << "creation of group " << _group << " failed";
		throw std::runtime_error(temp.str());
	}

	//---------------------------create dataspace---------------------------------
	_dataspace_id_data = H5Screate_simple(2, _dims, NULL);
	if (_dataspace_id_data < 0)
	{
		std::ostringstream temp;
		temp << "creation of dataspace for data failed";
		throw std::runtime_error(temp.str());
	}

	//-------------------------------------create property list interface------------
	_plist_id = H5Pcreate(H5P_DATASET_CREATE);
	if (_plist_id < 0)
	{
		std::ostringstream temp;
		temp << "creation of property list interface failed";
		throw std::runtime_error(temp.str());
	}
	//-----------------------------attributes in dense storage---------------------
	_status = H5Pset_attr_phase_change(_plist_id,0,0);
	if (_status < 0)
	{
		std::ostringstream temp;
		temp << "H5Pset_attr_phase_change failed";
		throw std::runtime_error(temp.str());
	}
	//---------------------------setup compression------------------------------
	_status = H5Pset_scaleoffset (_plist_id, H5Z_SO_INT, H5Z_SO_INT_MINBITS_DEFAULT);
	if (_status < 0)
	{
		std::ostringstream temp;
		temp << "H5Pset_scaleoffset failed";
		throw std::runtime_error(temp.str());
	}
	//-----------------------------setup compression chunks----------------------
	_status = H5Pset_chunk(_plist_id, 2, _cdims);
	if (_status < 0)
	{
		std::ostringstream temp;
		temp << "creation of chunks failed";
		throw std::runtime_error(temp.str());
	}

	//-----------------------------create dataset----------------------------------

	_dataset_id = H5Dcreate2(_file_id, _loc.c_str(), H5T_NATIVE_SHORT, _dataspace_id_data, H5P_DEFAULT, _plist_id, H5P_DEFAULT);
	if (_dataset_id < 0)
	{
		std::ostringstream temp;
		temp << "creation of dataset " << _loc << " failed";
		throw std::runtime_error(temp.str());
	}

	//----------------------create attributes------------------------------------

	_initialiseAttributeMap();
	if(_attributes_enabled){
		std::string temp_loc=_group+"/attributes";
		_group_attr_id = H5Gcreate2(_file_id, temp_loc.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (_group_attr_id < 0)
		{
			std::ostringstream temp;
			temp << "creation of group " << temp_loc << " failed";
			throw std::runtime_error(temp.str());
		}
	}
	for (auto it : _attributes) {
		_adims = getDim(it.second);
		_dataspace_id_attr = H5Screate_simple(1, &_adims, NULL);
		if (_dataspace_id_attr < 0) {
			throw std::runtime_error("could not create dataspace while creating attributes");
		}

		hid_t plist_id_attr = H5Pcreate(H5P_DATASET_CREATE);
		if (plist_id_attr < 0)
		{
			std::ostringstream temp;
			temp << "creation of property list interface failed while creating attributes";
			throw std::runtime_error(temp.str());
		}
		//---------------------------setup compression------------------------------
		/**_status = H5Pset_scaleoffset (plist_id_attr, H5Z_SO_INT, H5Z_SO_INT_MINBITS_DEFAULT);
		if (_status < 0)
		{
			std::ostringstream temp;
			temp << "H5Pset_scaleoffset failed while creating attributes";
			throw std::runtime_error(temp.str());
		}
		//-----------------------------setup compression chunks----------------------
		_status = H5Pset_chunk(plist_id_attr, 1, &_adims);
		if (_status < 0)
		{
			std::ostringstream temp;
			temp << "creation of chunks failed while creating attributes";
			throw std::runtime_error(temp.str());
		}
		**/
		std::string temp_loc=_group+"/attributes/"+it.first;
		hid_t temp = H5Dcreate2(_file_id, temp_loc.c_str(), getType(it.second), _dataspace_id_attr, H5P_DEFAULT,plist_id_attr, H5P_DEFAULT);
		//hid_t temp = H5Acreate(_dataset_id, it.first.c_str(), getType(it.second), _dataspace_id_attr, H5P_DEFAULT, H5P_DEFAULT);
		if (temp < 0) {
			std::ostringstream tempError;
			tempError<<"could not create attribute "<<it.first;
			throw std::runtime_error(tempError.str());
		}
		setAttribute(it.second, temp);
		_attributes [it.first] = it.second;
		H5Pclose(plist_id_attr);
		H5Sclose(_dataspace_id_attr);
	}
}
	const std::string _filename;
	bool _open = false;
	bool _exists = false;
	bool _attributes_enabled=true;
	bool _minimal_attributes=false;
	std::string _plane = "horizontal";
	std::string _beam = "B1";
	std::string _group = "/" + _beam;
	std::string _loc = _group + "/" + _plane;
	std::size_t _turns = 4096;
	std::size_t _bunches = 3564;

	herr_t _status;

	//HDF file identifiers
	hsize_t _dims[2] = { _turns, _bunches };
	hid_t _file_access_plist_id = 0;
	hid_t _file_id = 0;
	hid_t _group_id = 0;
	hid_t _group_attr_id = 0;
	hid_t _root_group_id = -1;
	hid_t _dataspace_id_data = 0;
	hid_t _plist_id = 0;
	hid_t _dataset_id = 0;
	hsize_t _cdims[2] = { _turns, 1 };

	//attributes
	hid_t _dataspace_id_attr = 0;
	hsize_t _adims = 1;

	std::map<std::string, Attribute> _attributes;

	//diagnostics
	std::size_t _time_to_write_data;



};

}