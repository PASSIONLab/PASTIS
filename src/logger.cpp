/**
 * @file
 *	logger.cpp
 *
 * @author
 *	Oguz Selvitopi
 *
 * @date
 *
 * @brief
 *	Simple logging utility
 *
 * @todo
 *
 * @note
 *	
 */

#include <memory>

#include "../inc/logger.hpp"



namespace
pastis
{

std::shared_ptr<Logger>	 Logger::logger_ = nullptr;

}
