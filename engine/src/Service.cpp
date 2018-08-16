#include "Service.h"

#include "RotatorService.h"
#include "CameraObserver.h"

using namespace opa::math::game;

OPA_NAMESPACE_DECL2(opa, engine)

Service::Service() { m_game = 0; }

void Service::init(Game *game) {
    opa::utils::Initable::init();
    m_game = game;
}

ServiceManager::ServiceManager() { m_game = 0; }

void ServiceManager::fini() {
    for (auto &x : m_services)
        if (x.ND.enabled)
            x.ND.service->stop();
    m_services.clear();
}

void ServiceManager::init(Game *game) {
    m_game = game;
    register_service(RotatorService::KEY, ServicePtr(new RotatorService));
    register_service(CameraObserver::KEY, ServicePtr(new CameraObserver));
}

void ServiceManager::enable_service(const ServiceKey &key) {
    OPA_CHECK0(m_services.count(key));
    auto &status = m_services[key];
    OPA_CHECK0(status.enabled == false);
    status.enabled = true;
    status.service->start();
}

void ServiceManager::disable_service(const ServiceKey &key) {
    OPA_CHECK0(m_services.count(key));
    auto &status = m_services[key];
    OPA_CHECK0(status.enabled == true);
    status.service->stop();
    status.enabled = false;
}

void ServiceManager::register_service(const ServiceKey &key,
                                      ServicePtr service) {
    OPA_CHECK0(m_services.count(key) == 0);
    service->init(m_game);
    ServiceStatus status;
    status.enabled = false;
    status.service = service;
    m_services[key] = status;
}

OPA_NAMESPACE_DECL2_END
