import React from 'react';
import clsx from 'clsx';

function Feature({ Svg, title, description }) {
  return (
    <div className={clsx('col col--4')}>
      {Svg && (
        <div className="text--center">
          <Svg className={styles.featureSvg} role="img" />
        </div>
      )}
      <div className={clsx('text--center', 'padding-horiz--md')}>
        {title && <h3>{title}</h3>}
        {description && <p>{description}</p>}
      </div>
    </div>
  );
}

export default Feature;
