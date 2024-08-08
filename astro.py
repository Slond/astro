import sys
from datetime import datetime, timezone

import astropy.units as u
import ephem
import pytz

from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
from geopy.geocoders import Nominatim
from timezonefinder import TimezoneFinder


planets = ['Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
geolocator = Nominatim(user_agent='myapplication')

try:
    num_stars = int(sys.argv[3])
except IndexError:
    num_stars = 25

def get_planet_coordinates(planet_name, coords, time):
    observer = ephem.Observer()
    observer.lat = coords[0]
    observer.lon = coords[1]
    observer.elev = 0
    observer.date = time

    if planet_name == 'Mercury':
        planet = ephem.Mercury(observer)
    elif planet_name == 'Venus':
        planet = ephem.Venus(observer)
    elif planet_name == 'Mars':
        planet = ephem.Mars(observer)
    elif planet_name == 'Jupiter':
        planet = ephem.Jupiter(observer)
    elif planet_name == 'Saturn':
        planet = ephem.Saturn(observer)
    elif planet_name == 'Uranus':
        planet = ephem.Uranus(observer)
    elif planet_name == 'Neptune':
        planet = ephem.Neptune(observer)
    else:
        return None

    return {
        'name': planet.name,
        'ra': float(planet.ra),
        'dec': float(planet.dec),
        'alt': float(planet.alt),
        'az': float(planet.az)
    }


def star_list(num_stars: int, data, coords, time):
    location = EarthLocation(lat=coords[0] * u.deg, lon=coords[1] * u.deg, height=0 * u.m)
    
    if len(data) == 0:
        return "Данные о звездах не найдены."
    else:
        stars = data[0]
        stars.sort('Vmag')

        custom_simbad = Simbad()
        custom_simbad.add_votable_fields('flux(V)')

        count = 0
        for star in stars:
            if count >= num_stars:
                break
            ra = star['_RAJ2000'] * u.deg
            dec = star['_DEJ2000'] * u.deg
            coord = SkyCoord(ra=ra, dec=dec, frame='icrs')

            # Преобразование в горизонтальную систему координат (AltAz)
            altaz = coord.transform_to(AltAz(obstime=time, location=location))

            # Угол возвышения > 0 (над горизонтом)
            if altaz.alt > 0 * u.deg:
                result_table = custom_simbad.query_region(coord)
                if result_table:
                    name = result_table[0]['MAIN_ID']
                else:
                    name = "Unknown"
                    
                # print(f"Звезда: {name}, RA={ra}, Dec={dec}, Угловая величина: Alt={altaz.alt:.2f}, Az={altaz.az:.2f}, Видимая величина: {star['Vmag']:.2f}")
                print(f"Звезда: {name}, Угловая величина: Alt={altaz.alt:.2f}, Az={altaz.az:.2f}, Видимая величина: {star['Vmag']:.2f}")
                count += 1


def planet_list(coords, time):
    for planet in planets:
        planet_data = get_planet_coordinates(planet, coords, time)
        if planet_data and planet_data['alt'] > 0:
            # print(f"Планета: {planet_data['name']}, RA={planet_data['ra']}, Dec={planet_data['dec']}, Угловая величина: Alt={planet_data['alt']:.2f}, Az={planet_data['az']:.2f}")
            print(f"Планета: {planet_data['name']}, Угловая величина: Alt={planet_data['alt']:.2f}, Az={planet_data['az']:.2f}")


def main(city: str, time: datetime, Vmag_num: float, stars_num: int):
    location = geolocator.geocode(city, language='ru-RU')
    coords = location[1]

    observ_time = time  # Год, месяц, день, час, минута, секунда
    time = Time(observ_time)

    vizier = Vizier(columns=['_RAJ2000', '_DEJ2000', 'Vmag'])
    result = vizier.query_constraints(catalog='I/239/hip_main', Vmag='<'+str(Vmag_num))

    star_list(stars_num, result, coords, time)
        
    print('')
    
    planet_list(coords, observ_time)
    

def first_data(city, time):
    location = geolocator.geocode(city, language='ru-RU')
    coords = location[1]
    
    tf = TimezoneFinder()
    corrected_time = datetime.now(pytz.timezone(tf.timezone_at(lat=coords[0], lng=coords[1])))
    
    print(location)
    print(corrected_time.strftime("%d %B %Y %I:%M%p"))
    print()


if __name__ == '__main__':
    # time = datetime(2024, 8, 8, 11, 30, 0) # Можно указать время вручную
    time = datetime.now(timezone.utc)  # Текущее время
    city = sys.argv[1]
    first_data(city, time)
    main(city, time, sys.argv[2], num_stars)  # Город, время, порог видимой звездной величины, максимальное количество выводимых звезд
