#pragma once

#include <opa/engine/conf.h>
#include <opa/engine/Game.h>

OPA_NAMESPACE_DECL2(opa, engine)

class Service : public opa::utils::Initable {
  public:
    virtual void init(Game *game);
    Service();

    virtual void start() = 0;
    virtual void stop() = 0;

    OPA_ACCESSOR_PTR(Camera, m_game->camera(), camera)
    OPA_ACCESSOR_PTR(Scene, m_game->scene(), scene)
    OPA_ACCESSOR_PTR(InputHandler, m_game->input_handler(), input_handler)
    OPA_ACCESSOR_PTR(Game, m_game, game);

  protected:
    virtual void init() override {}
    Game *m_game;
};
OPA_DECL_SPTR(Service, ServicePtr);

typedef std::string ServiceKey;

class ServiceManager : public opa::utils::Initable {
  public:
    struct ServiceStatus {
        ServicePtr service;
        bool enabled;
    };
    ServiceManager();
    virtual void init(Game *game);
    virtual void fini() override;

    void enable_service(const ServiceKey &key);
    void disable_service(const ServiceKey &key);

  protected:
    Game *m_game;

  private:
    virtual void init() override {}

    void register_service(const ServiceKey &key, ServicePtr service);
    std::map<ServiceKey, ServiceStatus> m_services;
};
OPA_DECL_SPTR(ServiceManager, ServiceManagerPtr);

template <class TData, class TSub> class PubSub : public Updater {
  public:
    virtual void init(Game *game, UpdateMark mark) {
        m_game = game;
        m_game->updater().add(mark, UpdaterPtr(this));
    }
    virtual void update() {
        // will update sub
        for (auto &x : m_subs)
            x->notify(m_data);
    }
    virtual void push(TData data) { m_data = data; }

    void add_sub(std::shared_ptr<TSub> sub) { m_subs.insert(sub); }

  private:
    Game *m_game;
    std::set<std::shared_ptr<TSub> > m_subs;
    TData m_data;
};

OPA_NAMESPACE_DECL2_END
