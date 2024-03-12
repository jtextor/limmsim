#ifndef __ENTITY_H

#define __ENTITY_H

class Entity{
  public:
    virtual bool operator== (const Entity &) { return true; }
    virtual int hashValue() = 0;

    virtual Entity * clone() = 0;
};


#endif
